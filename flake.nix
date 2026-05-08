{
  description = "Neuro 120 final project";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    pyproject-nix = {
      url = "github:pyproject-nix/pyproject.nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs =
    { nixpkgs, pyproject-nix, ... }:
    let
      inherit (nixpkgs) lib;
      forAllSystems = lib.genAttrs lib.systems.flakeExposed;

      project = pyproject-nix.lib.project.loadPyproject {
        projectRoot = ./.;
      };

      pythonAttr = "python312";
    in
    {
      devShells = forAllSystems (system: {
        default =
          let
            pkgs = nixpkgs.legacyPackages.${system};
            python = pkgs.${pythonAttr};
            pythonEnv = python.withPackages (project.renderers.withPackages { inherit python; });
          in
          pkgs.mkShell { packages = [ pythonEnv ]; };
      });

      packages = forAllSystems (
        system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
          python = pkgs.${pythonAttr};
        in
        {
          default = python.pkgs.buildPythonPackage (project.renderers.buildPythonPackage { inherit python; });
        }
      );
    };
}
