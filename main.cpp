#include "headers2d.hpp"

int main(int argc, char **argv){

  if(argc < 2){
    cout << "No setup file given.\nExiting program.\n";
    exit(1);
  }

  setupAide setup( argv[1] );

  return 0;
}
