#include "Ising/SquareIsingModel.h"
#include "MonteCarlo.hpp"
#include "Spin2d/SquareXYModel.h"
#include "SimulatorDisplay/Animator.h"
#include <chrono>


std::shared_ptr<MonteCarloSimulator> prepare_model() {
  Params params;
  params.emplace("system_size", 256);
  params.emplace("J", 1.0);
  params.emplace("B", 0.0);
  params.emplace("Bp", 0.0);
  params.emplace("temperature", 1.0);
  params.emplace("cluster_update", 0);

  return std::make_shared<SquareXYModel>(params, 1);
}

#ifdef EMSCRIPTEN
void main_loop(void* arg) {
  Animator* animator = static_cast<Animator*>(arg);
  if (!animator->is_paused()) {
    animator->simulator->timesteps(1);
  }
  FrameData frame_data = animator->new_frame(animator->simulator->get_texture());
  for (auto key : frame_data.keys) {
    animator->simulator->key_callback(key);
  }

  if (frame_data.status_code == 0) {

  }
}

int main(int argc, char** argv) {
  Animator animator(prepare_model());
  int fps = 60;
  emscripten_set_main_loop_arg(main_loop, &animator, fps, 1);
  return 0;
}

#else
int main_loop(void* arg) {
  auto start = std::chrono::high_resolution_clock::now();

  Animator* animator = static_cast<Animator*>(arg);
  if (!animator->is_paused()) {
    animator->simulator->timesteps(1);
  }
  FrameData frame_data = animator->new_frame(animator->simulator->get_texture());
  for (auto key : frame_data.keys) {
    animator->simulator->key_callback(key);
  }

  auto end = std::chrono::high_resolution_clock::now();
  double dt = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count(); 
  int fps = 60;
  double target_dt = 1.0/fps; 
  if (dt < target_dt) { 
    sleep(target_dt - dt); 
  }

  return frame_data.status_code;
}

int main(int argc, char** argv) {
  Animator animator(prepare_model());
  animator.start(900, 900);
  int status = 1;
  while (status == 1) {
    status = main_loop(&animator);
  }

  return 0;
}

#endif
