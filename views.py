from django.shortcuts import render_to_response, get_object_or_404
from django.http import *

def index(request):
    return render_to_response('dismod_welcome.html', {})
