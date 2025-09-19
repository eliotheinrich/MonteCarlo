#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Shader {
  public:
    unsigned int ID;

    Shader()=default;

    Shader(const char* source, const char* vertex_path, const char* fragment_path);

    void use() const;

    void set_bool(const std::string &name, bool value) const;
    void set_int(const std::string &name, int value) const;
    void set_float(const std::string &name, float value) const;
    void set_float3(const std::string &name, float v1, float v2, float v3) const;
    void set_float4(const std::string &name, float v1, float v2, float v3, float v4) const;

  private:
    void check_compile_errors(unsigned int shader, std::string type);
};
