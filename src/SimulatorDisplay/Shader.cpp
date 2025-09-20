#include "Shader.h"

#include <format>
#include <GLES3/gl3.h>

Shader Shader::make() {
#ifdef EMSCRIPTEN
  const char* vertex_shader_code = R"(#version 300 es
precision mediump float;
in vec3 aPos;
in vec2 aTexCoord;
out vec2 TexCoord;
void main() {
    gl_Position = vec4(aPos, 1.0);
    TexCoord = aTexCoord;
})";

const char* fragment_shader_code = R"(#version 300 es
precision mediump float;
in vec2 TexCoord;
out vec4 FragColor;
uniform sampler2D tex;
void main() {
    FragColor = texture(tex, TexCoord);
})";
#else
  const char* vertex_shader_code = "#version 330 core\nlayout (location = 0) in vec3 aPos;\nlayout (location = 1) in vec2 aTexCoord;\nout vec2 TexCoord;\nvoid main() { gl_Position = vec4(aPos, 1.0); TexCoord = aTexCoord; };";
  const char* fragment_shader_code = "#version 330 core\nout vec4 FragColor;\nin vec2 TexCoord;\nuniform sampler2D tex;\nvoid main() {\nFragColor = texture(tex, TexCoord);\n}";
#endif

  unsigned int vertex, fragment;
  Shader shader;

  vertex = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertex, 1, &vertex_shader_code, NULL);
  glCompileShader(vertex);
  shader.check_compile_errors(vertex, "VERTEX");

  fragment = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragment, 1, &fragment_shader_code, NULL);
  glCompileShader(fragment);
  shader.check_compile_errors(fragment, "FRAGMENT");

  shader.ID = glCreateProgram();
  glAttachShader(shader.ID, vertex);
  glAttachShader(shader.ID, fragment);
  glLinkProgram(shader.ID);
  shader.check_compile_errors(shader.ID, "PROGRAM");

  glDeleteShader(vertex);
  glDeleteShader(fragment);
  return shader;
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
