from django.urls import path
from . import views

urlpatterns = [
    #homepage returns a react component that works as a SPA
    path("", views.home, name="home"),
    path("detail/<int:id>", views.home, name="home"),

    #below endpoints act as API endpoints and return JSON
    path("api/request", views.new_request, name="new_request"),
    path("api/requests", views.get_requests, name="get_requests"),
    path("api/status/<int:alignment_request_id>", views.status, name="status"),
    path("api/detail/<int:alignment_request_id>", views.alignment_detail, name="alignment_detail")

]
