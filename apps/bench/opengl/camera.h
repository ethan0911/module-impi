#pragma once
#ifndef OSPRAY_CAMERA_H
#define OSPRAY_CAMERA_H
#include "common.h"
#include "trackball.h"

class Camera {
private:
  vec2f mouse2screen(int x, int y, float width, float height) {
    return vec2f(2.0f * (float) x / width - 1.0f, 2.0f * (float) y / height - 1.0f);
  }

private:
  size_t width = 2048, height = 2048;
  float aspect = (float) width / height;
  float zNear = 1.f, zFar = 50.f;
  float fovy = 30.f;
  vec3f eye = vec3f(0.f, 0.f, -15.f); // this trackball requires camera to be
  vec3f focus = vec3f(0.f);           // initialized on negtive z axis with
  vec3f up = vec3f(0.f, 1.f, 0.f);   // y axis as the initial up vector !!!!
  Trackball ball;
  // OSPRay
  OSPCamera ospCamera = nullptr;
public:

  void Clean() {
    if (ospCamera != nullptr) {
      ospRelease(ospCamera);
      ospCamera = nullptr;
    }
  }

  void Init(OSPCamera camera, float scale) {
    if (camera == nullptr) { throw std::runtime_error("empty camera found"); }
    eye.z *= scale;
    ospCamera = camera;
    CameraUpdateView();
    CameraUpdateProj(this->width, this->height);
  }

  OSPCamera OSPRayPtr() { return this->ospCamera; }

  size_t CameraWidth() { return this->width; }

  size_t CameraHeight() { return this->height; }

  float CameraZNear() { return this->zNear; }

  float CameraZFar() { return this->zFar; }

  void CameraBeginZoom(float x, float y) {
    vec2f p = mouse2screen(x, y, this->width, this->height);
    this->ball.BeginZoom(p.x, p.y);
  }

  void CameraZoom(float x, float y) {
    vec2f p = mouse2screen(x, y, this->width, this->height);
    this->ball.Zoom(p.x, p.y);
    CameraUpdateView();
  }

  void CameraBeginDrag(float x, float y) {
    vec2f p = mouse2screen(x, y, this->width, this->height);
    this->ball.BeginDrag(p.x, p.y);
  }

  void CameraDrag(float x, float y) {
    vec2f p = mouse2screen(x, y, this->width, this->height);
    this->ball.Drag(p.x, p.y);
    CameraUpdateView();
  }

  void CameraUpdateView() {
    auto dir = -xfmVector(this->ball.Matrix().l, this->eye - this->focus);
    auto up  =  xfmVector(this->ball.Matrix().l, this->up);
    auto pos = (-dir + this->focus);
    ospSetVec3f(ospCamera, "pos", (osp::vec3f &) pos);
    ospSetVec3f(ospCamera, "dir", (osp::vec3f &) dir);
    ospSetVec3f(ospCamera, "up", (osp::vec3f &) up);
    ospCommit(ospCamera);
  }

  void CameraUpdateProj(size_t width, size_t height) {
    this->aspect = width / (float) height;
    this->width = width;
    this->height = height;
    ospSetf(ospCamera, "aspect", this->aspect);
    ospSetf(ospCamera, "fovy", this->fovy);
    ospCommit(ospCamera);
  }
};

#endif //OSPRAY_CAMERA_H
