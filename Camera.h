#ifndef CAMERA_H
#define CAMERA_H

#include <glad/glad.h>  
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>  // Added so the header knows what glm::mat4 is

// Tell the linker that these variables exist in Camera.cpp
extern glm::mat4 projection;
extern glm::mat4 view;

void mouse_callback(GLFWwindow* window, double xposIn, double yposIn);
void updateCameraMatrices(int displayWidth, int displayHeight);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

#endif