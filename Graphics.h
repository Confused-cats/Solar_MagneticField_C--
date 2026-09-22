#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <vector>

// Forward declarations/Includes
#include "Camera.h" 

// Shaders (GLSL sources & Shader compilation)
extern const char* vertexShaderSource;
extern const char* fragmentShaderSource;
extern const char* flatVertexShaderSource;
extern const char* flatFragmentShaderSource;

extern const char* computeShaderSource;
GLuint createComputeShaderProgram(const char* csSource);

extern bool isStreaming;
GLuint createShaderProgram(const char* vsSource, const char* fsSource);

// Mesh Wrapper Classes
struct Grid {
    GLuint VAO = 0, VBO = 0;
    int vertexCount = 0;

    GLint projLoc = -1, viewLoc = -1, modelLoc = -1, colorLoc = -1;

    void setup(float size, int divisions);
    void draw(GLuint shader, const glm::mat4& proj, const glm::mat4& viewMat);
    void cleanup();
};

struct StarSphere {
    GLuint VAO = 0, VBO = 0, EBO = 0;
    int indexCount = 0;

    GLint projLoc = -1, viewLoc = -1, modelLoc = -1, colorLoc = -1;

    void setup(float radius, int slices, int stacks);
    void draw(GLuint shader, const glm::mat4& proj, const glm::mat4& viewMat);
    void cleanup();
};

// Window & Input Management
GLFWwindow* StartGLFW();
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);

// Visual Utilities (Interacts with Physics) 
void generateFieldLines(std::vector<float>& lineVertices);