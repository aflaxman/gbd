""" Script to extract all data associated with a given disease model
to a json file, so that it can be moved to a new edition of dismod."""

from django.core import serializers
from dismod_data_server.models import *

dm_id = 2764

def move_model(dm_id=2764):
    dm = DiseaseModel.objects.get(id=dm_id)

    models = [ dm ]
    models += dm.params.all()
    models += dm.data.all()


    data = serializers.serialize("json", models)
    out = open("dm_%d.json"%dm_id, "w")
    out.write(data)
    out.close()

    print """now use loaddata to load the json file::
        python manage.py loaddata dm_*.json"""
