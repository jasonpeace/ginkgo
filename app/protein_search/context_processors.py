from django.conf import settings


def reactjs_assets_paths(request):
    staticfiles_base = settings.STATICFILES_BASE
    build_files = settings.REACT_JS_BUILD_DIR
    js_paths = [str(x.relative_to(staticfiles_base)) for x in build_files.glob("*.js")]
    css_paths = [str(x.relative_to(staticfiles_base)) for x in build_files.glob("*.css")] 
    return {
        "reactjs_assets_js_paths":js_paths,
        "reactjs_assets_css_paths":css_paths,
    }

