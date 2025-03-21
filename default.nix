{
  packageName,
  packageVersion,

  buildPythonPackage,

  # build-system
  gcc,
  setuptools_scm,
  cython,
}:

let
in buildPythonPackage rec {
  pname = packageName;
  version = packageVersion;
  pyproject = true;
  build-system = [ setuptools_scm ];
  src = ./.;

  propagatedBuildInputs = [
    gcc
    cython
  ];

  # pythonImportsCheck = [ packageName ];  # we currently build using python and not sage, thus import checking is not possible

  meta = {
    description = "SAGE/Python functions useful for studying S-boxes and Boolean functions such as computing the DDT, computing the Walsh spectrum, affine equivalence testing...";
    homepage = "https://github.com/lpp-crypto/sboxU";
  };
}

