
{
  description = "Development flake for Neuro 120 project";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    pyproject-nix = {
      url = "github:pyproject-nix/pyproject.nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs =
    { nixpkgs, pyproject-nix, ... }:
    let
      project = pyproject-nix.lib.project.loadPyproject {
        projectRoot = ./.;
      };
      pkgs = nixpkgs.legacyPackages.${system};
      python = pkgs.python312;
      system = "aarch64-darwin";
      localPkg = python.pkgs.buildPythonPackage (
        project.renderers.buildPythonPackage { inherit python; }
      );
    in
    {
      devShells.${system}.default =
        let
          pythonEnv = python.withPackages (
            ps: project.renderers.withPackages { inherit python; } ps
          );
        in
        pkgs.mkShell {
          packages = [ pythonEnv ];
          shellHook = ''
            export PYTHONPATH="$PWD:$PYTHONPATH"
          '';
        };
      packages.${system}.default = localPkg;
    };
}

