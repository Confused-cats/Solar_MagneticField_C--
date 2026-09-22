#include "Graphics.h"
#include "Physics.h"
#include "Camera.h"
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>


// Shader Source definition

const char* vertexShaderSource = R"glsl(
    #version 450 core
    layout (location = 0) in vec3 aPos;   
    layout (location = 1) in vec3 aColor; 

    out vec3 ParticleColor;

    uniform mat4 projection;
    uniform mat4 view;

    void main() {
        ParticleColor = aColor;
        gl_Position = projection * view * vec4(aPos, 1.0);
    }
)glsl";

const char* fragmentShaderSource = R"glsl(
    #version 450 core
    in vec3 ParticleColor;
    out vec4 FragColor;

    void main() {
        float dist = gl_FragCoord.z / gl_FragCoord.w;
        float fog = exp(-0.01 * dist) * 1.5f;
        FragColor = vec4(ParticleColor * fog, 1.0);
    }
)glsl";

const char* flatVertexShaderSource = R"glsl(
    #version 450 core
    layout (location = 0) in vec3 aPos;

    uniform mat4 projection;
    uniform mat4 view;
    uniform mat4 model;

    void main() {
        gl_Position = projection * view * model * vec4(aPos, 1.0);
    }
)glsl";

const char* flatFragmentShaderSource = R"glsl(
    #version 450 core
    out vec4 FragColor;
    uniform vec3 uColor;

    void main() {
        FragColor = vec4(uColor, 1.0);
    }
)glsl";

const char* computeShaderSource = R"glsl(
    #version 450 core

    layout(local_size_x = 256, local_size_y = 1, local_size_z = 1) in;

    struct ParticleState {
        vec4 pos;    // xyz = physics pos (m), w = active flag
        vec4 vel;    // xyz = velocity (m/s),  w = q/m ratio
        vec4 color;  // rgb = color, a = particle type ID
    };

    struct OutputVertex {
        vec4 pos;    
        vec4 color;  
    };

    layout(std430, binding = 0) buffer ParticleBuffer {
        ParticleState particles[];
    };

    layout(std430, binding = 1) buffer OutputVertexBuffer {
        OutputVertex lineVertices[];
    };

    uniform float uPhysicsToGraphics;

    void main() {
        uint gid = gl_GlobalInvocationID.x;
        if (gid >= particles.length()) return;

        uint lineIdx = gid * 2;
        ParticleState p = particles[gid];

        if (p.pos.w <= 0.0) {
            lineVertices[lineIdx].pos = vec4(0.0);
            lineVertices[lineIdx].color = vec4(0.0);
            lineVertices[lineIdx + 1].pos = vec4(0.0);
            lineVertices[lineIdx + 1].color = vec4(0.0);
            return;
        }

        vec3 pos = p.pos.xyz;
        vec3 vel = p.vel.xyz;

        // Transform real-world physics coordinates (meters) to graphics space
        vec3 gfxPos = pos * uPhysicsToGraphics;
        float speed = length(vel);
        vec3 velDir = speed > 0.0 ? vel / speed : vec3(0.0);
        vec3 gfxTail = gfxPos - (velDir * 0.05);

        // Head Vertex
        lineVertices[lineIdx].pos = vec4(gfxPos, 1.0);
        lineVertices[lineIdx].color = p.color;

        // Tail Vertex
        lineVertices[lineIdx + 1].pos = vec4(gfxTail, 1.0);
        lineVertices[lineIdx + 1].color = vec4(p.color.rgb * 0.12, 1.0);
    }
)glsl";

// Rendering Helper Implementations

GLuint createShaderProgram(const char* vsSource, const char* fsSource) {
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vsSource, NULL);
    glCompileShader(vertexShader);

    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fsSource, NULL);
    glCompileShader(fragmentShader);

    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    return shaderProgram;
}

GLuint createComputeShaderProgram(const char* csSource) {
    GLuint computeShader = glCreateShader(GL_COMPUTE_SHADER);
    glShaderSource(computeShader, 1, &csSource, NULL);
    glCompileShader(computeShader);

    // Error logging check
    GLint success;
    glGetShaderiv(computeShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetShaderInfoLog(computeShader, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::COMPUTE::COMPILATION_FAILED\n" << infoLog << std::endl;
    }

    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, computeShader);
    glLinkProgram(shaderProgram);
    glDeleteShader(computeShader);

    return shaderProgram;
}
//Grid Visual Wrapper 

void Grid::setup(float size, int divisions) {
    std::vector<float> vertices;
    float step = (2.0f * size) / divisions;

    // Scaling factor converting graphics space coordinates to meters
    const double GRAPHICS_TO_PHYSICS = STAR_RADIUS / 1.25;
    
    // Visual tuning factor: Adjust this to make the gravity funnel deeper or shallower
    const double POTENTIAL_SCALE = 1.0e-11; 

    // Helper lambda to calculate downward displacement based on Gravitational Potential
    auto getWarpedY = [&](float x, float z) -> float {
        double physX = static_cast<double>(x) * GRAPHICS_TO_PHYSICS;
        double physZ = static_cast<double>(z) * GRAPHICS_TO_PHYSICS;
        double physR = std::sqrt(physX * physX + physZ * physZ);

        double baseHeight = -1.5; // Base height offset for the grid
        double dip = 0.0;

        if (physR >= STAR_RADIUS) {
            // Gravitational Potential = GM / r
            dip = (G * STAR_MASS / physR) * POTENTIAL_SCALE;
        } else {
            // Cap the dip at the star's surface so it doesn't form an infinite spike| and idk man, gibberish dsg lyigvsykgkuyfk
            dip = (G * STAR_MASS / STAR_RADIUS) * POTENTIAL_SCALE;
        }

        return static_cast<float>(baseHeight - dip);
    };

    // Generate grid lines running parallel to the Z axis
    for (int i = 0; i <= divisions; ++i) {
        float x = -size + i * step;
        for (int j = 0; j < divisions; ++j) {
            float z1 = -size + j * step;
            float z2 = -size + (j + 1) * step;

            float y1 = getWarpedY(x, z1);
            float y2 = getWarpedY(x, z2);

            vertices.push_back(x);  vertices.push_back(y1); vertices.push_back(z1);
            vertices.push_back(x);  vertices.push_back(y2); vertices.push_back(z2);
        }
    }

    // Generate grid lines running parallel to the X axis
    for (int j = 0; j <= divisions; ++j) {
        float z = -size + j * step;
        for (int i = 0; i < divisions; ++i) {
            float x1 = -size + i * step;
            float x2 = -size + (i + 1) * step;

            float y1 = getWarpedY(x1, z);
            float y2 = getWarpedY(x2, z);

            vertices.push_back(x1); vertices.push_back(y1); vertices.push_back(z);
            vertices.push_back(x2); vertices.push_back(y2); vertices.push_back(z);
        }
    }

    vertexCount = static_cast<int>(vertices.size() / 3);

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);
}

void Grid::draw(GLuint shader, const glm::mat4& proj, const glm::mat4& viewMat) {
    glUseProgram(shader);
    
    if (projLoc == -1) {
        projLoc  = glGetUniformLocation(shader, "projection");
        viewLoc  = glGetUniformLocation(shader, "view");
        modelLoc = glGetUniformLocation(shader, "model");
        colorLoc = glGetUniformLocation(shader, "uColor");
    }

    glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(proj));
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(viewMat));
    
    glm::mat4 model = glm::mat4(1.0f);
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

    glUniform3f(colorLoc, 0.15f, 0.20f, 0.25f);

    glLineWidth(1.0f);
    glBindVertexArray(VAO);
    glDrawArrays(GL_LINES, 0, vertexCount);
    glBindVertexArray(0);
}

void Grid::cleanup() {
    if (VAO) glDeleteVertexArrays(1, &VAO);
    if (VBO) glDeleteBuffers(1, &VBO);
}

// Star Sphere Mesh

void StarSphere::setup(float radius, int slices, int stacks) {
    std::vector<float> vertices;
    std::vector<unsigned int> indices;
    const float PI = 3.1415926535f;

    for (int i = 0; i <= stacks; ++i) {
        float phi = PI * static_cast<float>(i) / stacks;
        for (int j = 0; j <= slices; ++j) {
            float theta = 2.0f * PI * static_cast<float>(j) / slices;

            float x = radius * sinf(phi) * cosf(theta);
            float y = radius * cosf(phi);
            float z = radius * sinf(phi) * sinf(theta);

            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);
        }
    }

    for (int i = 0; i < stacks; ++i) {
        for (int j = 0; j < slices; ++j) {
            unsigned int first = (i * (slices + 1)) + j;
            unsigned int second = first + slices + 1;

            indices.push_back(first);
            indices.push_back(second);
            indices.push_back(first + 1);

            indices.push_back(second);
            indices.push_back(second + 1);
            indices.push_back(first + 1);
        }
    }

    indexCount = static_cast<int>(indices.size());

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);
}

void StarSphere::draw(GLuint shader, const glm::mat4& proj, const glm::mat4& viewMat) {
    glUseProgram(shader); 

    if (projLoc == -1) {
        projLoc  = glGetUniformLocation(shader, "projection");
        viewLoc  = glGetUniformLocation(shader, "view");
        modelLoc = glGetUniformLocation(shader, "model");
        colorLoc = glGetUniformLocation(shader, "uColor");
    }

    glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(proj)); 
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(viewMat)); 

    glm::mat4 model = glm::rotate(glm::mat4(1.0f), static_cast<float>(starRotationAngle), glm::vec3(0.0f, 1.0f, 0.0f));
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model)); 

    glUniform3f(colorLoc, 0.8f, 0.45f, 0.15f); 

    glBindVertexArray(VAO); 
    glDrawElements(GL_TRIANGLES, indexCount, GL_UNSIGNED_INT, 0); 
    glBindVertexArray(0); 
}

void StarSphere::cleanup() {
    if (VAO) glDeleteVertexArrays(1, &VAO);
    if (VBO) glDeleteBuffers(1, &VBO);
    if (EBO) glDeleteBuffers(1, &EBO);
}

//Dynamic Field Lines

void generateFieldLines(std::vector<float>& lineVertices) {
    lineVertices.clear();

    const float starRadius = 1.25f; 
    const float stepSize = 0.04f;  
    const int numSteps = 800;      
    const int numLines = 20;      
    
    glm::vec3 lineColor = glm::vec3(0.0f, 0.6f, 1.0f);
    float rotationAngle = static_cast<float>(starRotationAngle);

    for (int i = 0; i < numLines; ++i) {
        float angle = (i * 2.0f * 3.1415926f) / numLines + rotationAngle;
        
        float latitudeOffset = 0.0f;
        if (i % 3 == 1) latitudeOffset = 0.12f; 
        if (i % 3 == 2) latitudeOffset = -0.12f;

        glm::vec3 seedPos = glm::vec3(
            cos(angle) * (starRadius + 0.05f),
            latitudeOffset,
            sin(angle) * (starRadius + 0.05f)
        );

        const double GRAPHICS_TO_PHYSICS = STAR_RADIUS / 1.25;

        // Along the field vector (+B)
        glm::vec3 currentPos = seedPos;
        for (int step = 0; step < numSteps; ++step) { 
            
            // Scale graphics coordinate UP to astronomical scale
            glm::dvec3 physicsPos = glm::dvec3(currentPos) * GRAPHICS_TO_PHYSICS;
            
            //Calculate field using the real physics coordinates
            glm::vec3 B = glm::vec3(calculateMagneticField(physicsPos));
            
            float bLength = glm::length(B);
            
            if (bLength < 1e-8f) break;

            glm::vec3 dir = B / bLength;
            glm::vec3 nextPos = currentPos + dir * stepSize;

            if (glm::length(nextPos) < starRadius - 0.02f || glm::length(nextPos) > 32.0f) { 
                break; 
            }

            lineVertices.push_back(currentPos.x); lineVertices.push_back(currentPos.y); lineVertices.push_back(currentPos.z);
            lineVertices.push_back(lineColor.r);  lineVertices.push_back(lineColor.g);  lineVertices.push_back(lineColor.b);
            lineVertices.push_back(nextPos.x);    lineVertices.push_back(nextPos.y);    lineVertices.push_back(nextPos.z);
            lineVertices.push_back(lineColor.r);  lineVertices.push_back(lineColor.g);  lineVertices.push_back(lineColor.b);

            currentPos = nextPos;
        }

        // Against the field vector (-B)
        currentPos = seedPos;
        for (int step = 0; step < numSteps; ++step) { 
            
            // 1. Scale up for the negative trace as well
            glm::dvec3 physicsPos = glm::dvec3(currentPos) * GRAPHICS_TO_PHYSICS;
            glm::vec3 B = glm::vec3(calculateMagneticField(physicsPos));
            
            float bLength = glm::length(B);
            
            //threshold 
            if (bLength < 1e-8f) break;

            glm::vec3 dir = -B / bLength;
            glm::vec3 nextPos = currentPos + dir * stepSize;

            if (glm::length(nextPos) < starRadius - 0.02f || glm::length(nextPos) > 32.0f) { 
                break;
            }

            lineVertices.push_back(currentPos.x); lineVertices.push_back(currentPos.y); lineVertices.push_back(currentPos.z);
            lineVertices.push_back(lineColor.r);  lineVertices.push_back(lineColor.g);  lineVertices.push_back(lineColor.b);
            lineVertices.push_back(nextPos.x);    lineVertices.push_back(nextPos.y);    lineVertices.push_back(nextPos.z);
            lineVertices.push_back(lineColor.r);  lineVertices.push_back(lineColor.g);  lineVertices.push_back(lineColor.b);

            currentPos = nextPos;

        }
    }
}

// Window & Input Callbacks 

bool isStreaming = false;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
    if (key == GLFW_KEY_F && action == GLFW_PRESS) {
        isStreaming = !isStreaming; 
    }
    if (key == GLFW_KEY_P && action == GLFW_PRESS) {
        useParkerSpiral = !useParkerSpiral;
        std::cout << "Field Architecture Toggled! Parker Spiral: " 
                  << (useParkerSpiral ? "ENABLED" : "DISABLED (Pure Dipole mode)") << std::endl;
    }
    if (key == GLFW_KEY_E && action == GLFW_PRESS) {
        enableElectricField = !enableElectricField;
        std::cout << "Induced Electric Field Toggled! Co-rotation E-Field: " 
                  << (enableElectricField ? "ENABLED" : "DISABLED") << std::endl;
    }
}

