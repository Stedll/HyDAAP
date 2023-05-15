import numpy as np

def euler_to_matrix(alpha, beta, gamma):
  # Conversion from Euler angles in ZYZ convention to a rotation matrix
  R_z1 = np.array([[np.cos(alpha), -np.sin(alpha), 0],
                   [np.sin(alpha), np.cos(alpha), 0],
                   [0, 0, 1]])
  R_y = np.array([[np.cos(beta), 0, np.sin(beta)],
                  [0, 1, 0],
                  [-np.sin(beta), 0, np.cos(beta)]])
  R_z2 = np.array([[np.cos(gamma), -np.sin(gamma), 0],
                   [np.sin(gamma), np.cos(gamma), 0],
                   [0, 0, 1]])
  R = np.dot(R_z1, np.dot(R_y, R_z2))
  return R

def euler_to_matrix(array):
  #array[2]=array[2]-np.pi/2
  # Conversion from Euler angles in ZYZ convention to a rotation matrix
  R_z1 = np.array([[np.cos(array[0]), -np.sin(array[0]), 0],
                   [np.sin(array[0]), np.cos(array[0]), 0],
                   [0, 0, 1]])
  R_y = np.array([[np.cos(array[1]), 0, np.sin(array[1])],
                  [0, 1, 0],
                  [-np.sin(array[1]), 0, np.cos(array[1])]])
  R_z2 = np.array([[np.cos(array[2]), -np.sin(array[2]), 0],
                   [np.sin(array[2]), np.cos(array[2]), 0],
                   [0, 0, 1]])
  R = np.dot(R_z1, np.dot(R_y, R_z2))
  return R

def matrix_to_euler(R):
  # Conversion from a rotation matrix to Euler angles in ZYZ convention
  beta = np.arccos(R[2, 2])
  alpha = np.arctan2(R[1, 2], R[0, 2])
  gamma = np.arctan2(R[2, 1], -R[2, 0])
  return alpha, beta, gamma

