from django.shortcuts import render, redirect
from django.core.files.storage import FileSystemStorage

from .forms import FileForm

def set_nightmode(request):
    request.session['nightmode'] = request.session.get('nightmode', True)
    mode = request.session['nightmode']
    if request.method == 'POST':
        request.POST['nightmode']
        request.session['nightmode'] = not mode
    return

def start(request):
    try:
        set_nightmode(request)
    except:
        pass
    return render(request, 'geco/start.html', {})

def input_custom(request):
    try:
        set_nightmode(request)
    except:
        if request.method == "POST":
            if request.POST['input_choice'] == "file":
                form = FileForm(request.POST)
                if form.is_valid():
                    try:
                        inputfile = request.FILES['input_file']
                        fs = FileSystemStorage()
                        filename = fs.save(inputfile.name, inputfile)
                        uploaded_file_url = fs.url(filename)
                        url = ",".join(uploaded_file_url.strip().split("/"))
                        return redirect('file_context',
                                        uploaded_url=url)
                    except:
                        pass
    context = {
        'file_form' : FileForm()
    }
    return render(request, 'geco/input_custom.html', context)

def file_context(request, uploaded_url):
    try:
        set_nightmode(request)
    except:
        pass
    return render(request, 'geco/context.html', {
                                            'uploaded_url' : uploaded_url,
                                                 })
