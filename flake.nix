{
  description = "Neuro 120 final project";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-25.11-darwin";
  };

  outputs = inputs @ { self, ... }:
    let
      system = "aarch64-darwin";
      pkgs = import inputs.nixpkgs { inherit system; };
      python = pkgs.python313;
    in
    {
      devShells.${system}.default = pkgs.mkShellNoCC {
        packages = [
          (python.withPackages (ps: with ps; [
            pip
            jupyter
            jupyterlab-vim
            ipywidgets

            matplotlib
            numpy
            scipy
            pandas
          ]))
        ];
      };
    };
}
