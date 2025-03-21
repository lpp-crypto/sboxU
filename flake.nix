{
  description = "SAGE/Python functions useful for studying S-boxes and Boolean functions such as computing the DDT, computing the Walsh spectrum, affine equivalence testing...";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-24.05";
  };

  outputs = { self, nixpkgs }: let
    system = "x86_64-linux";
    pkgs = import nixpkgs { inherit system; };
    packageName = "sboxU";
    packageVersion = "1.3.0";
  in {

    packages."${system}" = {
      "${packageName}" = pkgs.python3Packages.callPackage ./default.nix {
        inherit packageName packageVersion;
      };

      default = self.packages."${system}"."${packageName}";

      overlay = (final: prev: {
        pythonPackagesExtensions = prev.pythonPackagesExtensions ++ [
          (final: prev: {
            "${packageName}" = self.packages."${system}".default;
          })
        ];
      });
    };

    devShell."${system}" = let
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ self.packages."${system}".overlay ];
        };

	customSage = pkgs.sage.override {
          requireSageTests = false;
          extraPythonPackages = ps: with ps; [
            ps."${packageName}"
	  ];
        };
        customPython = with pkgs; (pkgs.sage.withPackages (ps: [
          ps."${packageName}"
        ]));

        packages = [ customSage ];
    in pkgs.mkShell {
      inherit packages;
    };

  };
}

