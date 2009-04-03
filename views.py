from django.shortcuts import render_to_response, get_object_or_404
from django.http import *

def index(request):
    return HttpResponseRedirect('public/')
