import os
import pytest
import autostreamtree
import autostreamtree.functions as funcs


@pytest.fixture
def shapefile_path():
    base_path = os.path.dirname(autostreamtree.__file__)
    file_path = os.path.join(base_path, 'data', 'test.shp')
    return file_path


@pytest.fixture
def coords_path():
    base_path = os.path.dirname(autostreamtree.__file__)
    file_path = os.path.join(base_path, 'data', 'test.pointCoords.txt')
    return file_path
