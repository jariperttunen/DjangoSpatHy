import os
import subprocess
import glob
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
# Create your views here.
import SpathyDemo.spath as spath
from SpathyDemo.forms import SpathyFiguresForm

def index(request):
    return render(request,spath.mainpage,{'maintitle':spath.maintitle})

def run(request):
    spathy=spath.spathy_runtime
    sp=spath.spathy_dir
    subprocess.call([spath.spathy_runtime],cwd=spath.spathy_dir)
    figure_ls = glob.glob(spath.spathy_dir+'*.png')
    figure_ls = [os.path.basename(f) for f in figure_ls]
    return render(request,spath.resultpage,{'maintitle':spath.maintitle,
                                            'results':spath.results,
                                            'fig_ls':figure_ls})
