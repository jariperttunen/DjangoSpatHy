# Generated by Django 2.0.5 on 2018-05-15 13:33

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('SpathyDemo', '0003_auto_20180515_1332'),
    ]

    operations = [
        migrations.AlterField(
            model_name='spathyfiguresmodel',
            name='figure_1',
            field=models.FileField(default='name', max_length=255, upload_to=''),
        ),
        migrations.AlterField(
            model_name='spathyfiguresmodel',
            name='figure_2',
            field=models.FileField(default='name', max_length=255, upload_to=''),
        ),
    ]
