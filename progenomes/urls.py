from django.urls import path
from . import views
from . import api


urlpatterns = [
    # EGGNOG
    path('input/', views.input, name='input'),
    path('context/<str:query>/',
         views.context,
         name='context'),
    path('api/getcontext/<str:query>/',
         api.get_context,
         name='get_context'),
    path('api/tree/<str:query>/',
         api.get_tree,
         name='get_tree'),
    path('api/egglevels/',
         api.get_eggNOG_levels,
         name='get_eggNOG_levels'),
]
