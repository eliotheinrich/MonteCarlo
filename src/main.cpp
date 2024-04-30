#include "SquareIsingModel.h"
#include "SquareXYModel.h"
#include <Frame.h>

using namespace dataframe;

int main() {
  Params params;
  params.emplace("system_size", 16.0);
  params.emplace("J", 1.0);
  params.emplace("B", 0.0);
  params.emplace("Bp", 0.0);
  params.emplace("temperature", 1.0);
	

  //SquareIsingModel model(params, 1);
  //model.timesteps(100);


  SquareXYModel model2(params, 1);
  model2.timesteps(100);

}
