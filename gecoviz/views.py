from django.shortcuts import render, redirect
from django.core.files.storage import FileSystemStorage

from .forms import FileForm, NewickForm

def set_nightmode(request):
    request.session['nightmode'] = request.session.get('nightmode', False)
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
    return render(request, 'gecoviz/start.html', {})

def input_custom(request):
    def getURL(fs, inputFile):
        filename = fs.save(inputFile.name, inputFile)
        uploaded_file_url = fs.url(filename)
        url = ",".join(uploaded_file_url.strip().split("/"))
        return url
    try:
        set_nightmode(request)
    except:
        if request.method == "POST":
            if request.POST['input_choice'] == "file":
                form = FileForm(request.POST)
                if form.is_valid():
                    # try:
                    fs = FileSystemStorage()
                    inputFile = request.FILES['input_file']
                    file_url = getURL(fs, inputFile)
                    try:
                        inputNewick = request.FILES['newick_file']
                        newick_url = getURL(fs, inputNewick)
                    except:
                        newick_url = '_'
                    return redirect('file_context',
                                    file_url=file_url,
                                    newick_url=newick_url)
                    # except:
                        # pass
    context = {
        'newick_form' : NewickForm(),
        'file_form' : FileForm(),
    }
    return render(request, 'gecoviz/input_custom.html', context)

def file_context(request, file_url , newick_url):
    try:
        set_nightmode(request)
    except:
        pass
    return render(request, 'gecoviz/context.html', {
                                            'file_url' : file_url,
                                            'newick_url' : newick_url,
                                                 })

def documentation(request):
    try:
        set_nightmode(request)
    except:
        pass
    return render(request, 'gecoviz/documentation.html', {})
