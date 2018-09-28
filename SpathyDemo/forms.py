from django import forms
from SpathyDemo.models import SpathyFiguresModel

class SpathyFiguresForm(forms.ModelForm):
    class Meta:
        model = SpathyFiguresModel
        fields = ('figure_1','figure_2')
