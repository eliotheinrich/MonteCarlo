#include "MonteCarlo.hpp"
#include <memory>
#include <set>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <format>
#include <emscripten/emscripten.h>
#include <emscripten/html5.h>
#include <GLES3/gl3.h>


const char* vertex_shader_src = R"(
attribute vec2 a_position;
attribute vec2 a_texcoord;
varying vec2 v_texcoord;
void main() {
    v_texcoord = a_texcoord;
    gl_Position = vec4(a_position, 0.0, 1.0);
}
)";

const char* fragment_shader_src = R"(
precision mediump float;
varying vec2 v_texcoord;
uniform sampler2D u_texture;
void main() {
    gl_FragColor = texture2D(u_texture, v_texcoord);
}
)";

GLuint vbo = 0;
GLuint texture = 0;

class Shader {
  public:
    GLuint id;

    Shader(const char* src, GLenum type) {
      id = glCreateShader(type);
      glShaderSource(id, 1, &src, nullptr);
      glCompileShader(id);
      GLint ok = 0;
      glGetShaderiv(id, GL_COMPILE_STATUS, &ok);

      if (!ok) {
        GLint len = 0;
        glGetShaderiv(id, GL_INFO_LOG_LENGTH, &len);
        std::vector<char> log(len);
        glGetShaderInfoLog(id, len, nullptr, log.data());
        glDeleteShader(id);
        throw std::runtime_error(std::format("Shader error: {}", log.data()));
      }
    }
};

EM_BOOL key_callback(int eventType, const EmscriptenKeyboardEvent *e, void *userData);

class glProgram {
  private:
    GLuint id;

  public:
    std::set<std::string> keys_pressed;
    std::map<std::string, double> key_timer;

    Shader vs;
    Shader fs;
    std::shared_ptr<MonteCarloSimulator> simulator;
    int width;
    int height;

    int get_id() const {
      return id;
    }

    int get_width() const {
      return width;
    }

    int get_height() const {
      return height;
    } 

    void press(const std::string& key) {
      if (key == "Escape") {
        std::cout << "Registered escape!\n";
      } else if (key == " ") {
        std::cout << "Registered pause!\n";
      } else {
        simulator->key_callback(key);
      }
    }

    glProgram(const char* vs_src, const char* fs_src, std::shared_ptr<MonteCarloSimulator> simulator, int width, int height)
    : vs(vs_src, GL_VERTEX_SHADER), fs(fs_src, GL_FRAGMENT_SHADER), simulator(simulator), width(width), height(height) {
      id = glCreateProgram();
      glAttachShader(id, vs.id);
      glAttachShader(id, fs.id);
      glDeleteShader(vs.id);
      glDeleteShader(fs.id);
      glBindAttribLocation(id, 0, "a_position");
      glBindAttribLocation(id, 1, "a_texcoord");
      glLinkProgram(id);
      GLint ok = 0;
      glGetProgramiv(id, GL_LINK_STATUS, &ok);
      if (!ok) {
        GLint len = 0;
        glGetProgramiv(id, GL_INFO_LOG_LENGTH, &len);
        std::vector<char> log(len);
        glGetProgramInfoLog(id, len, nullptr, log.data());
        glDeleteProgram(id);
        throw std::runtime_error(std::format("Program link error: {}", log.data()));
      }

      emscripten_set_keydown_callback(EMSCRIPTEN_EVENT_TARGET_DOCUMENT, this, EM_TRUE, key_callback);
      emscripten_set_keyup_callback(EMSCRIPTEN_EVENT_TARGET_DOCUMENT, this, EM_TRUE, key_callback);
    }

    ~glProgram() {
      glDeleteProgram(id);
    }
};

EM_BOOL key_callback(int eventType, const EmscriptenKeyboardEvent *e, void *userData) {
  glProgram* program = static_cast<glProgram*>(userData);

  if (eventType == EMSCRIPTEN_EVENT_KEYDOWN) {
    program->keys_pressed.insert(e->key);
  } else if (eventType == EMSCRIPTEN_EVENT_KEYUP) {
    program->keys_pressed.erase(e->key);
  }

  return EM_TRUE;
}


glProgram init_gl(std::shared_ptr<MonteCarloSimulator> simulator) {
  // Initialize context and make current
  EmscriptenWebGLContextAttributes attr;
  emscripten_webgl_init_context_attributes(&attr);
  attr.majorVersion = 2; // Request WebGL2 for float textures
  EMSCRIPTEN_WEBGL_CONTEXT_HANDLE ctx = emscripten_webgl_create_context("#canvas", &attr);
  emscripten_webgl_enable_extension(ctx, "OES_texture_float");
  emscripten_webgl_enable_extension(ctx, "EXT_color_buffer_float");

  int w, h;
  emscripten_get_screen_size(&w, &h);
  EMSCRIPTEN_RESULT r = emscripten_set_canvas_element_size("#canvas", w, h);
  emscripten_webgl_make_context_current(ctx);

  glProgram program(vertex_shader_src, fragment_shader_src, simulator, w, h);

  // Fullscreen quad (two triangles), with texture coords
  const GLfloat vertices[] = {
    // x,   y,   u,   v
    -1.f,-1.f, 0.f, 0.f,
    1.f,-1.f, 1.f, 0.f,
    -1.f, 1.f, 0.f, 1.f,
    1.f, 1.f, 1.f, 1.f
  };

  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

  Texture tex = program.simulator->get_texture();

  // Generate texture from float RGBA array
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, tex.m, tex.n, 0,
      GL_RGBA, GL_FLOAT, tex.data());

  // Sampling settings
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  return program;
}

int frame = 0;

void draw_frame(void* arg) {
  glProgram* program = static_cast<glProgram*>(arg);
  glViewport(0, 0, program->get_width(), program->get_height());
  glClear(GL_COLOR_BUFFER_BIT);

  glUseProgram(program->get_id());
  glBindBuffer(GL_ARRAY_BUFFER, vbo);

  glEnableVertexAttribArray(0); // position
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (void*)0);

  glEnableVertexAttribArray(1); // texcoord
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (void*)(2 * sizeof(GLfloat)));
  
  program->simulator->timesteps(1);
  Texture tex = program->simulator->get_texture();
  for (auto k : program->keys_pressed) {
    program->press(k);
    std::cout << std::format("pressed k = {}\n", k);
  }
  
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, tex.n, tex.m,
      GL_RGBA, GL_FLOAT, tex.data());

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, texture);
  glUniform1i(glGetUniformLocation(program->get_id(), "u_texture"), 0);

  // Draw 2 triangles as a strip
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
}

#include "Ising/SquareIsingModel.h"
#include "Spin2d/SquareXYModel.h"

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

int main() {
  glProgram program = init_gl(prepare_model());
  int fps = 60;
  emscripten_set_main_loop_arg(draw_frame, &program, 0, EM_TRUE);
  return 0;
}

