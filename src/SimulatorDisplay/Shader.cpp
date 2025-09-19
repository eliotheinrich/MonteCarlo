#include "Shader.h"

#include <format>
#include <glad/glad.h>

Shader::Shader(const char* source, const char* vertex_path, const char* fragment_path) {
  std::string vertex_code;
  std::string fragment_code;
  std::ifstream vertex_shader_file;
  std::ifstream fragment_shader_file;

  vertex_shader_file.exceptions (std::ifstream::failbit | std::ifstream::badbit);
  fragment_shader_file.exceptions (std::ifstream::failbit | std::ifstream::badbit);
  try {
    std::string source_dir = std::format("{}/", source);
    vertex_shader_file.open(source_dir + vertex_path);
    fragment_shader_file.open(source_dir + fragment_path);
    std::stringstream vertex_shader_stream, fragment_shader_stream;

    vertex_shader_stream << vertex_shader_file.rdbuf();
    fragment_shader_stream << fragment_shader_file.rdbuf();

    vertex_shader_file.close();
    fragment_shader_file.close();

    vertex_code   = vertex_shader_stream.str();
    fragment_code = fragment_shader_stream.str();

  }

  catch (std::ifstream::failure& e) {
    std::cout << "ERROR::SHADER::FILE_NOT_SUCCESSFULLY_READ: " << e.what() << std::endl;
  }

  const char* vertex_shader_code = vertex_code.c_str();
  const char* fragment_shader_code = fragment_code.c_str();

  unsigned int vertex, fragment;

  vertex = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertex, 1, &vertex_shader_code, NULL);
  glCompileShader(vertex);
  check_compile_errors(vertex, "VERTEX");

  fragment = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragment, 1, &fragment_shader_code, NULL);
  glCompileShader(fragment);
  check_compile_errors(fragment, "FRAGMENT");

  ID = glCreateProgram();
  glAttachShader(ID, vertex);
  glAttachShader(ID, fragment);
  glLinkProgram(ID);
  check_compile_errors(ID, "PROGRAM");

  glDeleteShader(vertex);
  glDeleteShader(fragment);
}

void Shader::use() const { 
  glUseProgram(ID); 
}

void Shader::set_bool(const std::string &name, bool value) const {         
  glUniform1i(glGetUniformLocation(ID, name.c_str()), (int)value); 
}

void Shader::set_int(const std::string &name, int value) const { 
  glUniform1i(glGetUniformLocation(ID, name.c_str()), value); 
}

void Shader::set_float(const std::string &name, float value) const { 
  glUniform1f(glGetUniformLocation(ID, name.c_str()), value); 
}

void Shader::set_float3(const std::string &name, float v1, float v2, float v3) const {
  glUniform3f(glGetUniformLocation(ID, name.c_str()), v1, v2, v3); 
}

void Shader::set_float4(const std::string &name, float v1, float v2, float v3, float v4) const {
  glUniform4f(glGetUniformLocation(ID, name.c_str()), v1, v2, v3, v4); 
}

void Shader::check_compile_errors(unsigned int shader, std::string type) {
  int success;
  char infoLog[1024];
  if (type != "PROGRAM") {
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
      glGetShaderInfoLog(shader, 1024, NULL, infoLog);
      std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
    }
  }
  else {
    glGetProgramiv(shader, GL_LINK_STATUS, &success);
    if (!success) {
      glGetProgramInfoLog(shader, 1024, NULL, infoLog);
      std::cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
    }
  }
}
