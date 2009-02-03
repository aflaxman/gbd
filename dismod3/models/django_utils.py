def copy_model_instance(obj):
    """
    make a copy of a django model
    does not copy ManyToMany associations
    by miracle2k, http://www.djangosnippets.org/snippets/1040/
    """
    from django.db.models import AutoField

    initial = dict([(f.name, getattr(obj, f.name))
                    for f in obj._meta.fields
                    if not isinstance(f, AutoField) and\
                    not f in obj._meta.parents.values()])
    return obj.__class__(**initial)

