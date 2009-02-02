"""
utility functions that several views will use.
"""
import StringIO
import csv

import pylab as pl

MIMETYPE = {'png': 'image/png',
            'svg': 'image/svg+xml',
            'eps': 'application/postscript',
            'ps': 'application/postscript',
            'pdf': 'application/pdf',
            'csv': 'text/csv',
            'json': 'application/json',
            }

command_list = {'edit': ['edit'],
                'move': ['prev', 'next'],
                'sex': ['male', 'female', 'total'],
                'format': ['png', 'svg', 'pdf', 'csv', 'json'],
                }

id_delta = {'prev': -1, 'next': 1}

def template_params(object, **params):
    template_params = {'command_list': command_list, 'id': object.id, 'obj': object}
    template_params.update(params)
    return template_params

def figure_data(format):
    """
    return a string containing the representation of the current
    matplotlib figure.  format must be something matplotlib
    understands
    """
    f = StringIO.StringIO()
    pl.savefig(f, format=format)
    f.seek(0)
    return f.read()


def csv_str(headings, rows):
    """
    return a string containing a csv version of the table with the
    given headings and row data
    """
    f = StringIO.StringIO()
    csv_writer = csv.writer(f)
    csv_writer.writerow(headings)
    csv_writer.writerows(rows)
    f.seek(0)
    return f.read()

def clear_plot(width=4*1.5, height=3*1.5):
    fig = pl.figure(figsize=(width,height))
    pl.clf()
    return fig

def label_plot(title, **params):
    pl.xlabel('Age (years)', **params)
    pl.ylabel('Rate (per 1.0)', **params)
    pl.title(str(title), **params)

def objects_to_id_str(objs):
    return '_'.join([str(o.id) for o in objs])

def id_str_to_objects(id_str, obj_class):
    id_list = [int(id) for id in id_str.split('_')]
    return obj_class.objects.filter(id__in=id_list)
