#include <iostream>
#include <stdio.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm> 

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Camera.h"
#include "Physics.h"
#include "Graphics.h"



const double PHYSICS_TO_GRAPHICS = 1.25 / STAR_RADIUS;

GLFWwindow* StartGLFW();

int main() {
    srand(static_cast<unsigned>(time(0))); 

    GLFWwindow* window = StartGLFW();
    if (window == nullptr) return -1;

    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE); 

    GLuint particleShader = createShaderProgram(vertexShaderSource, fragmentShaderSource);
    GLuint flatShader = createShaderProgram(flatVertexShaderSource, flatFragmentShaderSource);


    Grid referenceGrid;
    referenceGrid.setup(18.0f, 60);

    StarSphere star;
    star.setup(1.25f, 60, 60);

    // Generate GPU buffers 
    GLuint VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    size_t maxVertices = 120000; 
    size_t bufferSize = maxVertices * 6 * sizeof(float);
    glBufferData(GL_ARRAY_BUFFER, bufferSize, nullptr, GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0); 

    // MAagnetic field lines buffer setup
    unsigned int fieldVAO, fieldVBO;
    glGenVertexArrays(1, &fieldVAO);
    glGenBuffers(1, &fieldVBO);

    glBindVertexArray(fieldVAO);
    glBindBuffer(GL_ARRAY_BUFFER, fieldVBO);

    // Dynamic allocation structure for trace line buffer
    glBufferData(GL_ARRAY_BUFFER, 16 * 800 * 2 * 12 * sizeof(float), NULL, GL_DYNAMIC_DRAW);

    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);

    // Setup callbacks
    glfwSetKeyCallback(window, key_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    float deltaTime = 0.0f; 
    float lastframe = 0.0f;
    static bool cPressedLastFrame = false;
    static bool vPressedLastFrame = false;

    const GLubyte* version = glGetString(GL_VERSION);
    std::cout << "OpenGL version " << version << std::endl;

    const GLubyte* glslversion = glGetString(GL_SHADING_LANGUAGE_VERSION);
    std::cout << "GLSL Version " << glslversion << "\n" << std::endl;


    std::cout << "CONTROLS" << std::endl;
    std::cout << "left Click and Drag to move the camera" << std::endl;
    std::cout << "P = Toggle Parker Spiral" << std::endl;
    std::cout << "E = Toggle Electric Field" << std::endl;
    std::cout << "C = Inject Cosmic Rays" << std::endl;
    std::cout << "V = Trigger Coronal Mass Ejection" << std::endl;
    std::cout << "B = Trigger Coronal Mass Ejection [CARRINGTON-CLASS]" << std::endl;

    GLint particleProjLoc = glGetUniformLocation(particleShader, "projection");
    GLint particleViewLoc = glGetUniformLocation(particleShader, "view");

    GLuint computeShader = createComputeShaderProgram(computeShaderSource);

    // Cache uniform locations outside the render loop
    GLint computeActiveCountLoc = glGetUniformLocation(computeShader, "uActiveCount");
    GLint computeScaleLoc       = glGetUniformLocation(computeShader, "uPhysicsToGraphics");

    const size_t MAX_PARTICLES = 100000;

    // Pre-allocate static/persistent transfer buffer to prevent runtime allocations
    std::vector<float> gpuParticleData(MAX_PARTICLES * 12);

    // Storage Buffer for particle physical state (q/m, position, velocity, color)
    GLuint particleSSBO;
    glCreateBuffers(1, &particleSSBO);
    glNamedBufferData(particleSSBO, MAX_PARTICLES * sizeof(float) * 12, nullptr, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, particleSSBO);

    // Storage Buffer for output line vertices 
    GLuint computeRenderVBO;
    glCreateBuffers(1, &computeRenderVBO);
    glNamedBufferData(computeRenderVBO, MAX_PARTICLES * 2 * 8 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, computeRenderVBO);

    // Configure computeVAO so glDrawArrays can render directly from computeRenderVBO
    GLuint computeVAO;
    glGenVertexArrays(1, &computeVAO);
    glBindVertexArray(computeVAO);
    glBindBuffer(GL_ARRAY_BUFFER, computeRenderVBO);

    // Position attribute (location = 0)
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Color attribute (location = 1)
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(4 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);

    while (!glfwWindowShouldClose(window)) {
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastframe;
        lastframe = currentFrame;

        if (deltaTime > 0.02f) deltaTime = 0.02f; 

        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS || isStreaming) {
            injectPlasmaField(plasmaParticles); 
        }

        if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS) {
            if (!cPressedLastFrame) {
                injectCosmicRay(plasmaParticles);
                cPressedLastFrame = true;
            }
        } else {
            cPressedLastFrame = false;
        }
        if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS) {
            if (!vPressedLastFrame) {
                std::cout <<  "\nCoronal Mass Ejection (CME) Ejecting\n" << std::endl;
                triggerCoronalMassEjection(plasmaParticles);
                vPressedLastFrame = true;
            }
        } else {
            vPressedLastFrame = false;
        }

        static bool bPressedLastFrame = false;

        if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS) {
            if (!bPressedLastFrame) {
                triggerCarringtonEvent(plasmaParticles);
                bPressedLastFrame = true;
            }
        } else {
            bPressedLastFrame = false;
        }

        updateParticles(plasmaParticles, deltaTime);

        GLuint numParticles = static_cast<GLuint>(plasmaParticles.size());

        if (numParticles > 0) {
            size_t activeCount = std::min(static_cast<size_t>(numParticles), MAX_PARTICLES);

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(activeCount); ++i) {
                const auto& p = plasmaParticles[i];
                size_t idx = static_cast<size_t>(i) * 12;

                gpuParticleData[idx + 0]  = static_cast<float>(p.position.x);
                gpuParticleData[idx + 1]  = static_cast<float>(p.position.y);
                gpuParticleData[idx + 2]  = static_cast<float>(p.position.z);
                gpuParticleData[idx + 3]  = p.active ? 1.0f : 0.0f;

                gpuParticleData[idx + 4]  = static_cast<float>(p.velocity.x);
                gpuParticleData[idx + 5]  = static_cast<float>(p.velocity.y);
                gpuParticleData[idx + 6]  = static_cast<float>(p.velocity.z);
                gpuParticleData[idx + 7]  = static_cast<float>(p.charge / p.mass);

                gpuParticleData[idx + 8]  = p.color.r;
                gpuParticleData[idx + 9]  = p.color.g;
                gpuParticleData[idx + 10] = p.color.b;
                gpuParticleData[idx + 11] = static_cast<float>(p.type);
            }

            glNamedBufferSubData(particleSSBO, 0, activeCount * 12 * sizeof(float), gpuParticleData.data());
        }

        glClearColor(0.02f, 0.02f, 0.03f, 1.0f); 
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        int width, height;
        glfwGetWindowSize(window, &width, &height);
        glViewport(0, 0, width, height);

        updateCameraMatrices(width, height); 
        
        referenceGrid.draw(flatShader, projection, view);
        star.draw(flatShader, projection, view);

        std::vector<float> fieldLineVertices;
        generateFieldLines(fieldLineVertices);

        if (!fieldLineVertices.empty()) {
            glUseProgram(particleShader);
            glUniformMatrix4fv(particleProjLoc, 1, GL_FALSE, glm::value_ptr(projection));
            glUniformMatrix4fv(particleViewLoc, 1, GL_FALSE, glm::value_ptr(view));

            glNamedBufferSubData(fieldVBO, 0, fieldLineVertices.size() * sizeof(float), fieldLineVertices.data());

            glBindVertexArray(fieldVAO);
            glLineWidth(2.2f);
            glDrawArrays(GL_LINES, 0, (GLsizei)(fieldLineVertices.size() / 6));
            glBindVertexArray(0);
        }

        // GPU Compute pass (Visualization / Vertex generation only)
        if (numParticles > 0) {
            glUseProgram(computeShader);
            glUniform1i(computeActiveCountLoc, static_cast<GLint>(numParticles));
            glUniform1f(computeScaleLoc, static_cast<float>(PHYSICS_TO_GRAPHICS));

            // Dispatch threads in workgroups of 256
            GLuint numWorkgroups = (numParticles + 255) / 256;
            glDispatchCompute(numWorkgroups, 1, 1);

            // Memory barrier: guarantee compute writes finish before drawing starts
            glMemoryBarrier(GL_VERTEX_ATTRIB_ARRAY_BARRIER_BIT | GL_SHADER_STORAGE_BARRIER_BIT);

            // GPU render pass (direct draw from ssbo [active particles only])
            glUseProgram(particleShader);
            glUniformMatrix4fv(particleProjLoc, 1, GL_FALSE, glm::value_ptr(projection));
            glUniformMatrix4fv(particleViewLoc, 1, GL_FALSE, glm::value_ptr(view));

            glBindVertexArray(computeVAO);
            glLineWidth(2.2f);
            glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(numParticles * 2));
            glBindVertexArray(0);
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    } 

    // Cleanup 
    referenceGrid.cleanup();
    star.cleanup();
    glDeleteProgram(flatShader);

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &fieldVAO);
    glDeleteBuffers(1, &fieldVBO);
    
    // Compute Shader & SSBOs
    glDeleteVertexArrays(1, &computeVAO);
    glDeleteBuffers(1, &particleSSBO);
    glDeleteBuffers(1, &computeRenderVBO);
    glDeleteProgram(computeShader);

    glDeleteProgram(particleShader);

    glfwTerminate();
    return 0;
}

GLFWwindow* StartGLFW() {
    if (!glfwInit()) {
        std::cerr << "Failed to start/initialize GLFW" << std::endl;
        return nullptr;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    #ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // Required for Apple ecosystems
    #endif

    GLFWwindow* window = glfwCreateWindow(800, 600, "3D Plasma Lorentz / Magnetic Field Sandbox", NULL, NULL);
    return window;
}