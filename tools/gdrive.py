"""
Functions for uploading data to the 
shared Google Drive analysis folder.
"""

from pathlib import Path
import gspread
from gspread_pandas import Spread, Client

# connect to Drive account for Gesenius Project
gc = gspread.service_account()
client = Client()

def get_existing_sheet(foldername, filename):
    """Retrieve existing Drive sheet if it exists.
    
    Args:
        foldername: str of the folder name
        filename: str of the filename
        
    Returns: 
        if sheet exists, dict of sheet data 
        with keys: id, name, path; otherwise 
        returns None
    """
    
    # NB: avoid double-naming folders in the drive
    # since it is not possible to provide a full path
    folders = client.find_folders(foldername)
    if not folders:
        return None
    
    # search for first matching sheet name
    # within the matched folders
    for folder in folders:
        for sheet in client.list_spreadsheet_files_in_folder(folder['id']):
            if sheet['name'] == filename:
                return sheet
        
    # no match; return a None
    return None

def get_sheet(foldername, filename, exist_ok=True):
    """Retrieves a sheet by foldername and file to perform operations on.
    
    Args:
        foldername: str name of the folder
        filename: str name of the sheet to be saved
    Returns:
        gcspread Spreadsheet object
    """
    
    # retrieve sheet if it already exists
    exist_sheet = get_existing_sheet(foldername, filename)
        
    # load sheet
    if exist_sheet:
        if not exist_ok:
            raise Exception('File already exists! Set exist_ok=True to force.')
        sheet = gc.open_by_key(exist_sheet['id']) 
    else:
        folders = client.find_folders(foldername)
        if not folders:
            raise Exception('Folder does not exist!')
        
        # NB: selects the first available folder with matching name
        # thus you must avoid double-naming folders!! There is no
        # way I know of r/n to do this with a full path.
        folder = folders[0]
        sheet = gc.create(filename, folder_id=folder['id'])
        
    return sheet

def df_to_drive(df, foldername, filename, worksheet_index=0, **spread_kwargs):
    """Write a DataFrame to Google Drive."""
    sh = get_sheet(foldername, filename)    
    spread = Spread(sh.id)
    skwargs = {'replace': True}
    skwargs.update(spread_kwargs)
    spread.df_to_sheet(df, **skwargs)

