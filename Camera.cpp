
// Camera angle viewing
#include <glad/glad.h>
#include "Camera.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <cmath>

glm::mat4 projection;
glm::mat4 view;

// camera state variables
static float yaw = 0.0f;
static float pitch = 0.0f;
static float camRadius = 20.0f;
static float lastX = 1200.0f;
static float lastY = 900.0f;
static bool firstMouse = true;

void mouse_callback(GLFWwindow* window, double xposIn, double yposIn) {
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
        
        if (firstMouse) {
            lastX = xpos;
            lastY = ypos;
            firstMouse = false;
        }

        float xoffset = xpos - lastX;
        float yoffset = lastY - ypos; 
        lastX = xpos;
        lastY = ypos;

        float sensitivity = 0.2f;
        xoffset *= sensitivity;
        yoffset *= sensitivity;

        yaw   -= xoffset;
        pitch -= yoffset;

        if (pitch > 89.0f)  pitch = 89.0f;
        if (pitch < -89.0f) pitch = -89.0f;
    } 
    else {
        firstMouse = true;
    }
}

void updateCameraMatrices(int displayWidth, int displayHeight) {
    if (displayHeight == 0) displayHeight = 1;
    float aspect = static_cast<float>(displayWidth) / displayHeight;

    // Handle Projection Matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    // Target global variable
    projection = glm::perspective(glm::radians(45.0f), aspect, 0.1f, 100.0f);
    glLoadMatrixf(glm::value_ptr(projection));

    
    glm::vec3 camPos;
    camPos.x = camRadius * cos(glm::radians(pitch)) * sin(glm::radians(yaw));
    camPos.y = camRadius * sin(glm::radians(pitch));
    camPos.z = camRadius * cos(glm::radians(pitch)) * cos(glm::radians(yaw));

    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    // Target the global variable
    view = glm::lookAt(camPos, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    glLoadMatrixf(glm::value_ptr(view));
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    camRadius -= static_cast<float>(yoffset) * 0.5f; 
    
    if (camRadius < 1.5f) camRadius = 1.5f; 
    if (camRadius > 50.0f) camRadius = 50.0f;
}