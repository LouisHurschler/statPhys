{
  description = "statPhys code";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs {
        inherit system;
      };
    in
    {
      formatter.${system} = pkgs.nixpkgs-fmt;

      devShell.${system} = pkgs.mkShell rec {
        nativeBuildInputs = with pkgs; [
        clang-tools
        tcsh
        ];
        buildInputs = with pkgs; [ 
        SDL2
        SDL2_image
        ];

        CPATH = pkgs.lib.makeSearchPathOutput "dev" "include" buildInputs;
        LD_LIBRARY_PATH = pkgs.lib.makeLibraryPath buildInputs;
      };
    };
}
