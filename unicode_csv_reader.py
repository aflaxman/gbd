import csv

def unicode_csv_reader(unicode_csv_data, dialect=csv.excel, **kwargs):
    """ csv.py doesn't do Unicode; encode temporarily as UTF-8: with this method."""
    csv_reader = csv.reader(utf_8_encoder(unicode_csv_data),
                            dialect=dialect, **kwargs)
    for row in csv_reader:
        # decode UTF-8 back to Unicode, cell by cell:
        yield [unicode(cell, 'utf-8') for cell in row]

def utf_8_encoder(unicode_csv_data):
    """ Helper method for unicode csv reader."""
    for line in unicode_csv_data:
        yield line.encode('utf-8')
