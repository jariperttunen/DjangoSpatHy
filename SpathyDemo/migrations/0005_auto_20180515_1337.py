# Generated by Django 2.0.5 on 2018-05-15 13:37

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('SpathyDemo', '0004_auto_20180515_1333'),
    ]

    operations = [
        migrations.AlterField(
            model_name='spathyfiguresmodel',
            name='figure_1',
            field=models.CharField(default='name', max_length=255),
        ),
        migrations.AlterField(
            model_name='spathyfiguresmodel',
            name='figure_2',
            field=models.CharField(default='name', max_length=255),
        ),
    ]
