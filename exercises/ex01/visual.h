//
// This file contains the variables and functions used
// for visualizing the trajectory of the two-dimensional
// ideal gas on the fly.
//
// This is cool, but not necessary to do the exercise.
// Also, it demands that you have installed the SDL
// graphics library.
//
// If the code is not compiled with the macro "VISUAL",
// these funcions are replaced by dummy functions that
// are called by the code, but do nothing.
//
#ifndef VISUAL
// "VISUAL" macro is not set - define dummy functions
void init_visual( ) { return; }
void draw_visual( double* ) { return; }
void kill_visual( ) { return; }
bool wait_visual( ) { return false; }
void wait_for_close( ) { return; }
#else
// "VISUAL" macro is set - we apply visualization using SDL
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL_timer.h>
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
// PARAMETERS OF THE WINDOW
//
const int window_edge      = 1000; // edge of the square window in points
const int window_margin    = 50;   // margin surrounding the square window
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
// VARIABLES AND PARAMETERS FOR THE VISUALIZATION
//
SDL_Window*   window;   // visualization window
SDL_Renderer* renderer; // visualization renderer
SDL_Surface*  surface;  // visualization object
SDL_Texture*  texture;  // visualization texture
SDL_Rect      position; // image position
SDL_Event     event;    // to detect a possible window close
//
const double box_to_window = window_edge / box_edge;   // system-to-window length conversion factor
const double rad_in_window = rad_atom * box_to_window; // atom radius mapped to its window size
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
void init_visual( ) {
  //
  // INITIALIZE THE VISUALIZATION
  //
  // returns zero on success else non-zero
  if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
    printf("error initializing SDL: %s\n", SDL_GetError());
  }
  //
  char txt[200];
  sprintf (txt, "2D IDEAL GAS - N=%d T=%.1f [K] B=%.1f [nm] A=%.1f [nm^2] - R=%.2f [nm] M=%.1f [g/mol] - Dt = %.2f [ps]",
	   num_atoms, ini_temp, box_edge, box_area, rad_atom, mas_atom, time_step);
  //
  // create a window
  window = SDL_CreateWindow(txt, 
			    SDL_WINDOWPOS_CENTERED,
			    SDL_WINDOWPOS_CENTERED,
			    window_edge + 2*window_margin, window_edge + 2*window_margin, 0);
  //
  // creates a renderer to render our images
  renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
  //
  // creates a surface to load an image into the main memory
  surface = IMG_Load("sphere.png");
  //
  // loads image to our graphics hardware memory.
  texture = SDL_CreateTextureFromSurface(renderer, surface);
  //
  // clears main-memory
  SDL_FreeSurface(surface);
  //
  // connects our texture with its position
  SDL_QueryTexture(texture, NULL, NULL, &position.w, &position.h);
  //
  // adjust height and width of our image box
  position.w = 2*rad_in_window;
  position.h = 2*rad_in_window;
  //
  return;
}
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
void draw_visual( double crd[] ) {
  //
  // DRAW THE ATOMS AT THEIR CURRENT POSITIONS AND THE WOX WALLS
  //
  // clear the rendering
  SDL_RenderClear(renderer);
  //
  // draw the atoms
  for ( int n = 0; n < num_atoms; n++ ) {
    position.x = window_margin + (window_edge/box_edge)*crd[2*n] - 0.5 * position.w;
    position.y = window_margin + (window_edge/box_edge)*crd[2*n+1] - 0.5 * position.h;
    SDL_RenderCopy(renderer, texture, NULL, &position);    
  }
  //  
  // draw the box walls
  SDL_SetRenderDrawColor(renderer, 0,0,0, SDL_ALPHA_OPAQUE);
  SDL_RenderDrawLine(renderer, window_margin, window_margin, window_margin, window_margin+window_edge);
  SDL_RenderDrawLine(renderer, window_margin, window_margin+window_edge, window_margin+window_edge, window_margin+window_edge);
  SDL_RenderDrawLine(renderer, window_margin+window_edge, window_margin+window_edge, window_margin+window_edge, window_margin);
  SDL_RenderDrawLine(renderer, window_margin+window_edge, window_margin, window_margin, window_margin);  
  SDL_SetRenderDrawColor(renderer, 255,255,255, SDL_ALPHA_OPAQUE);
  //
  // apply the the rendering
  SDL_RenderPresent(renderer);
  //
  return;
  //
}
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
bool wait_visual( ) {
  //
  // WAIT BEFORE MOVING TO NEXT FRAME AND CHECK IN CASE USER CLICKS OFF WINDOW
  //
  // calculates to 60 fps
  SDL_Delay(1000 / 60);
  //
  // events management
  SDL_PollEvent(&event);
  //
  // handling of close button
  if ( event.type == SDL_QUIT ) { return true; }
  //
  return false;
}
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
void kill_visual( ) {
  //
  // TERMINATE THE VISUALIZATION
  //
  // destroy texture
  SDL_DestroyTexture(texture);
  // destroy renderer
  SDL_DestroyRenderer(renderer);
  // destroy window
  SDL_DestroyWindow(window);
  // close SDL
  SDL_Quit();
  //
  return;
}
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
void wait_for_close( ) {
  //
  // WAIT THAT THE USER CLICKS OFF THE WINDOW
  //
  while ( !SDL_PollEvent(&event) || event.type != SDL_QUIT );  
  //
  return; 
}
//
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//
// endif for the macro "VISUAL" being defined
//
#endif
