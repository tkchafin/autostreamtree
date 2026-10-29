import importlib
import pkgutil


def test_import_all_modules():
    # Import all submodules in the resistnet package
    package_name = 'autostreamtree'
    package = importlib.import_module(package_name)

    for importer, modname, ispkg in pkgutil.walk_packages(
            package.__path__, package_name + '.'):
        importlib.import_module(modname)


def test_network_module_is_light():
    # autostreamtree.network (used by other packages, e.g. ResistNet) needs
    # only the network libraries: no VCF reading or Mantel tests
    import subprocess
    import sys
    code = ("import sys, autostreamtree.network; "
            "print(sorted({'pysam', 'mantel'} & set(sys.modules)))")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True,
                         text=True, check=True).stdout.strip()
    assert out == "[]"
