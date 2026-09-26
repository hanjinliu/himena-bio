import pytest
from himena.testing import install_plugin


@pytest.fixture(scope="session", autouse=True)
def init_pytest(request):
    install_plugin("himena-bio")
