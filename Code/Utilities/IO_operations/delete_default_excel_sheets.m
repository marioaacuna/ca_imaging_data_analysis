function delete_default_excel_sheets(filename, varargin)

p = inputParser();
addParameter(p, 'sheet_name', 'Sheet')  % EN: Sheet, DE: Tabelle, etc. (Lang. dependent)
parse(p, varargin{:})
sheet_name = p.Results.sheet_name;

% Open Excel file.
objExcel = actxserver('Excel.Application');
objExcel.Workbooks.Open(filename);
% Delete sheets
try
    % Throws an error if the sheets do not exist.
    objExcel.ActiveWorkbook.Worksheets.Item([sheet_name '1']).Delete;
    objExcel.ActiveWorkbook.Worksheets.Item([sheet_name '2']).Delete;
    objExcel.ActiveWorkbook.Worksheets.Item([sheet_name '3']).Delete;
end
% Save, close and clean up.
objExcel.ActiveWorkbook.Save;
objExcel.ActiveWorkbook.Close;
objExcel.Quit;
objExcel.delete;

