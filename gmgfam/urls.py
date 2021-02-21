from django.urls import path
from . import views
from . import api

urlpatterns = [
    # GMGFAM
    path('input/', views.input, name='input'),
    path('api/getcontext/<datatype>/<query>/<cutoff>/',
         api.get_context,
         name='get_context'),
    path('api/tree/<query>/',
         api.get_tree,
         name='get_tree'),
    path('api/egglevels/',
         api.get_eggNOG_levels,
         name='get_eggNOG_levels'),
]
