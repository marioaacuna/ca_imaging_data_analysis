import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    import h5py
import scipy.io as spio


def load(filename, force_dictionary=False, return_metadata=False, **kwargs):
    """
    this function should be called instead of direct spio.load
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """
    # Remove default parameters
    user_kwargs = kwargs.keys()
    if 'struct_as_record' in user_kwargs:
        kwargs.pop('struct_as_record')
    if 'squeeze_me' in user_kwargs:
        kwargs.pop('squeeze_me')

    try:
        data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True, **kwargs)
        # Convert objects to nested dictionaries and return
        data = _matobj_check_keys(data)
        # Extract metadata
        metadata = {v: data[v] for v in ['__version__', '__header__', '__globals__']}
        # Remove metadata
        data = {v: data[v] for v in data.keys() if v not in ['__version__', '__header__', '__globals__']}

    except NotImplementedError:
        with h5py.File(filename, mode='r') as f:
            # Read top-level
            data = dict({k: f[k] for k in list(f.keys()) if not k.startswith('#')})
            # Convert objects to nested dictionaries and return
            data = _hdf5_check_keys(data)
        # No metadata are available here
        metadata = dict()

    # If user requested only some variables, return only those
    variable_names = kwargs.get('variable_names', None)
    if variable_names is not None:
        # Make 'variable_names' a list
        if isinstance(variable_names, str):
            variable_names = [variable_names]
        data = dict({v: data[v] for v in variable_names})

    else:
        variable_names = list(data.keys())

    # Return result directly if only one variable requested
    if len(variable_names) == 1 and not force_dictionary:
        data = data[variable_names[0]]

    if not return_metadata:
        return data
    else:
        return data, metadata


def _matobj_check_keys(dict):
    """Checks if entries in dictionary are mat-objects. If yes, _matobj_to_dict is
    called to change them to nested dictionaries.
    """
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _matobj_to_dict(dict[key])

    return dict


def _matobj_to_dict(matobj):
    """A recursive function which constructs from matobjects nested dictionaries."""
    dictionary = dict()
    for fieldname in matobj._fieldnames:
        element = matobj.__dict__[fieldname]
        if isinstance(element, spio.matlab.mio5_params.mat_struct):
            dictionary[fieldname] = _matobj_to_dict(element)
        else:
            dictionary[fieldname] = element

    return dictionary


def _hdf5_check_keys(dictionary):
    """Checks if entries in dictionary are mat-objects. If yes, _matobj_to_dict is
    called to change them to nested dictionaries.
    """
    for key in dictionary:
        if isinstance(dictionary[key], h5py._hl.group.Group):
            dictionary[key] = _hdf5_to_dict(dictionary[key])

    return dictionary


def _hdf5_to_dict(hdf5_group):
    """A recursive function which constructs from hdf5 groups nested dictionaries."""
    dictionary = dict()
    for fieldname in hdf5_group.keys():
        element = hdf5_group[fieldname]
        if isinstance(element, h5py._hl.group.Group):
            dictionary[fieldname] = _hdf5_to_dict(element)
        else:
            dictionary[fieldname] = element[:]

    return dictionary
