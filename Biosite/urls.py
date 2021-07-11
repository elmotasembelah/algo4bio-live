from django.urls import path
from . import views


urlpatterns = [
    path('home/', views.homepage, name='home'),

    path('naiveExactMatching/', views.naiveExactMatchingPage, name='naive exact matching'),

    path('boyerMooreExactMatching/', views.boyer_MooreExactMatchingPage, name='boyer moore exact matching'),

    path('naiveApproximateMatching/', views.naiveApproximateMatchingPage, name='naive approximate matching'),

    path('boyerMooreApproximatematching/', views.boyerMooreApproximatematchingPage, name='boyer moore approximate matching'),

    path('hammingDistance/', views.hammingDistancePage, name='hamming distance'),

    path('editDistance/', views.editDistancePage, name='edit distance'),

    path('binarySearchKmerIndexing/', views.binarySearchKmerIndexingPage, name='binary search kmer indexing'),

    path('hashmapKmerIndexing/', views.hashmapKmerIndexingPage, name='hashmap kmer indexing'),

    path('localAlignment/', views.localAlignmentPage, name='localAlignment'),

    path('globalAlignment/', views.globalAlignmentPage, name='global alignment'),

]
