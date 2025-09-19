#include "Animator.h"
#include <unistd.h>
#include <map>

#include <functional>
#include <future>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#define SIMULATOR_DISPLAY_PATH "/Users/eliotheinrich/Projects/MonteCarlo/src/SimulatorDisplay"

static void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
  glViewport(0, 0, width, height);
}

class Animator::WindowImpl {
  public:
    GLFWwindow* impl;
    WindowImpl(GLFWwindow* impl) : impl(impl) {}
};

class Animator::InputProcessImpl {
  private:
    using KeyMap = std::map<int, std::function<void(Animator&, GLFWwindow*)>>;
    KeyMap keymap;
    std::map<int, double> prevtime;

    std::vector<int> key_buffer;

  public:
    InputProcessImpl() {
      keymap[GLFW_KEY_SPACE] = [](Animator& sim, GLFWwindow* window) { sim.paused = !sim.paused; };
      keymap[GLFW_KEY_ESCAPE] = [](Animator& sim, GLFWwindow* window) { glfwSetWindowShouldClose(window, true); };
    }

    double update_key_time(int key_code, double time) {
      if (!prevtime.contains(key_code)) {
        prevtime[key_code] = 0.0;
      }

      double t = time - prevtime[key_code];
      prevtime[key_code] = time;
      return t;
    }

    auto get_key_callback() {
      auto callback = [](GLFWwindow *window, int key, int scancode, int action, int mods) {
        auto& self = *static_cast<Animator*>(glfwGetWindowUserPointer(window));
        if (action == GLFW_RELEASE) {
          if (self.input_processor->keymap.contains(key)) {
            auto func = self.input_processor->keymap[key];
            func(self, window);
          }

          self.input_processor->key_buffer.push_back(key);
        }
      };

      return callback;
    }

    const std::vector<int> get_key_buffer() const {
      return key_buffer;
    }

    void clear_key_buffer() {
      key_buffer.clear();
    }
};

Animator::Animator(std::shared_ptr<MonteCarloSimulator> simulator, const Color& background_color) : simulator(simulator), background_color(background_color) {
  // Setting up state variables
  paused = false;

  frame_data = Texture();
  input_processor = std::make_unique<InputProcessImpl>();
}

Animator::~Animator() {
  glDeleteVertexArrays(1, &VAO);
  glDeleteBuffers(1, &VBO);
  glDeleteBuffers(1, &EBO);

  glfwTerminate();
}

void Animator::init_buffers() {
  // Setting up GL variables
  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);
  glGenBuffers(1, &EBO);
  glGenTextures(1, &texture_idx); 

  // Configure vertex attributes
  glBindVertexArray(VAO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);

  float vertices[] = {
    // positions         // frame_data coords
    1.0f,  1.0f, 0.0f,   1.0f, 1.0f,   // top right
    1.0f, -1.0f, 0.0f,   1.0f, 0.0f,   // bottom right
    -1.0f, -1.0f, 0.0f,   0.0f, 0.0f,   // bottom left
    -1.0f,  1.0f, 0.0f,   0.0f, 1.0f    // top left 
  };

  unsigned int indices[] = {
    0, 1, 3, // first triangle
    1, 2, 3  // second triangle
  };

  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_DYNAMIC_DRAW);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_DYNAMIC_DRAW);

  // Configure textures
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, texture_idx);  
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

  // Unbind buffers
  glBindVertexArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindTexture(GL_TEXTURE_2D, 0);  

  // Activate shader
  std::cout << "Loading shaders at path " << SIMULATOR_DISPLAY_PATH << "\n";
  shader = Shader(SIMULATOR_DISPLAY_PATH, "vertex_texture_shader.vs", "fragment_texture_shader.fs");
  shader.set_int("tex", 0);

  input_processor = std::make_unique<InputProcessImpl>();
}

void Animator::start(size_t width, size_t height) {
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
  glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_FALSE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

  GLFWwindow* glfw_window = glfwCreateWindow(width, height, "Simulator", NULL, NULL);
  window = std::make_unique<WindowImpl>(glfw_window);
  if (window->impl == NULL) {
    std::cout << "Failed to create GLFW window" << std::endl;
    glfwTerminate();
    return;
  }
  glfwMakeContextCurrent(window->impl); 

  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return;
  }  

  init_buffers();

  glViewport(0, 0, width, height);
  glfwSetFramebufferSizeCallback(window->impl, framebuffer_size_callback);  

  glfwSetWindowUserPointer(window->impl, this);
  auto key_callback = input_processor->get_key_callback();
  glfwSetKeyCallback(window->impl, key_callback);
}

FrameData Animator::new_frame(const Texture& new_data) {
  if (!paused) {
    frame_data = new_data;
    glClearColor(background_color.r, background_color.g, background_color.b, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    // Draw frame_data
    shader.use();
    //glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBindTexture(GL_TEXTURE_2D, texture_idx);  
    glActiveTexture(GL_TEXTURE0);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, frame_data.n, frame_data.m, 0, GL_RGBA, GL_FLOAT, (void*) frame_data.data());

    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glBindTexture(GL_TEXTURE_2D, 0);

    glfwSwapBuffers(window->impl);
  }
  
  glfwPollEvents();

  int status_code = glfwWindowShouldClose(window->impl) ? 0 : 1;
  std::vector<int> keys = input_processor->get_key_buffer();
  input_processor->clear_key_buffer();
  return {status_code, keys};
}

FrameData Animator::new_frame(const std::vector<float>& new_data, size_t n, size_t m) {
  Texture texture(new_data, n, m);
  return new_frame(texture);
}
