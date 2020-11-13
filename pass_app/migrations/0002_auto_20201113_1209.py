# Generated by Django 3.0.5 on 2020-11-13 12:09

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pass_app', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='compound',
            name='model_structure',
        ),
        migrations.AddField(
            model_name='compound',
            name='name',
            field=models.CharField(default='', max_length=200),
        ),
        migrations.AddField(
            model_name='compound',
            name='structure',
            field=models.CharField(default='', max_length=200),
        ),
    ]
