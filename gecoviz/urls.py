from django.urls import path
from . import api
from . import views

urlpatterns = [
    path('', views.start, name='start'),
    # path(r'api/info/<str:search_type>/<str:query>/',
         # api.info),
    # path(r'api/newick/<str:query>/', api.newick),
    # path(r'api/context/<str:query>/<int:cutoff>/', api.context),
    # CUSTOM INPUT
    path('custom/input/', views.input_custom, name='input_custom'),
    path('custom/filecontext/<str:uploaded_url>/',
         views.file_context,
         name='file_context'),
]

