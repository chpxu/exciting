{
  description = "Opinionated Flake for Python/C/C++/JS Development";
  inputs.nixpkgs = {
    url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  };
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = {
    self,
    nixpkgs,
    flake-utils,
  }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = import nixpkgs {
        inherit system;
        config.allowUnfree = true;
        config.allowUnfreePredicate = _: true;
        overlays = [
          (import ./nix/overlays)
        ];
      };
      # See nix/packages/default.nix to see what to pass into attrs set
      attrs = {
        inherit pkgs;
        lib = pkgs.lib;
        installC = false;
        installPython = true;
        installJS = false;
        installLatex = false;
        useLLVM = false;
        llvmVer = "17";
        pythonVer = "311";
        nodeVer = "20";
        enableVSCodeSetup = true;
      };
      # imports the shell and package configuration from `nix/default.nix`
      configuration = import ./nix attrs;
    in {
      # Set formatter
      formatter = pkgs.alejandra;
      devShells.default = configuration.shellOverride {
        nativeBuildInputs = [
          pkgs.bashInteractive
          pkgs.pkg-config
          (
            if attrs.useLLVM
            then pkgs."llvmPackages_${attrs.llvmVer}".libstdcxxClang
            else pkgs.gcc
          )
        ];
        packages = configuration.packages ++ [pkgs.bashInteractive];
      };
    });
}
