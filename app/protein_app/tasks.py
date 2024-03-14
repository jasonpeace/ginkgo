from celery import shared_task
from protein_app.modules.alignment import find_first_match


@shared_task()
def alignment_match(query: str, path: str, alignment_request_id: int):
    """
    This task will call the find_first_match function and is intended to be
    run by the celery instance.

    """
    find_first_match(query, path, alignment_request_id)
