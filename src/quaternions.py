"""Plate Rotations in Terms of Quaternions

Rotating plates using quaternions
"""

__author__ = 'Tobias Stephan'

import numpy as np
import math
import quaternion

def euler2quat(x):
    """Read data from R
    
    Returns the quaternion of an Euler pole and Euler angle
    
    """
    
    e = np.array([x['x'], x['y'], x['z']])  # / np.linalg.norm(np.array([x['x'], x['y'], x['z']]))
    q_sc = math.cos(x["angle"]/2) # scalar part
    q_vec = e * math.sin(x["angle"]/2) # vector part
    
    #q = q_sc + q_vec # unit real quaternion
    R = np.quaternion(q_sc, q_vec[0], q_vec[1], q_vec[2])
    
    return R

  
def euler_angle(R1, R2):
    """Euler angle from a concetanation of two rotations
    
    Returns Euler angle (in radians) associated with the rotation R2 following R1 
    
    
    Parameters
    ----------
    R1, R2 : quaternion or array of quaternions
             The quaternion(s) need not be normalized, but must be nonzero
    
    Returns
    -------
    w : float array
        Angle in radians
    
    """
    w = 2*math.acos((R2.w*R1.w) - np.dot(quaternion.as_vector_part(R2), quaternion.as_vector_part(R1)))
    
    return(w)

def euler_axis(R1, R2, w=None):
    """Euler axis from a concetanation of two rotations
    
    Returns Euler axis (in Cartesian coordinates) associated with the rotation R2 following R1 
    
    Parameters
    ----------
    R1, R2 : quaternion or array of quaternions
             The quaternion(s) need not be normalized, but must be nonzero. 
             R2 is the rotation that follows R1.
    w : None, float, or array of floats
        Euler angle
    
    Returns
    -------
    e : float array
        Euler axis
    
    """
  
    if w==None:
      w = euler_angle(R1, R2)
      
    e = 1 / math.sin(w/2) * (R2.w*quaternion.as_vector_part(R1) + R1.w*quaternion.as_vector_part(R2) + np.cross(quaternion.as_vector_part(R2), quaternion.as_vector_part(R1)))
    
    return e

def relative_rotation(R1, R2):
    """Relative rotation from two giving absolute rotations"""
    R = R2*R1.conjugate()
  
    return R


def rotate_vector_quat(u, q):
  """Rotation of vector
    
  Rotate vector u by quaternion q 
    
  Parameters
  ----------
  u : float array
      unit vector
    
  q : quaternion
    
  Returns
  -------
  w : float array
      rotated vector
    
  """
    
  #w = q*u*q.conjugate()
  m = quaternion.as_rotation_matrix(q)
  
  w = np.dot(m, u.T).T 
  
  return w
  


#def check_plate_circuit(R1, R2)


# Absolute rotation or point(s) P by rotation R
# Pprime = R*P

# Relative rotation of P1 wrt P2 with P1 rotating around R1 and P2 rotating around R2
# R1 * R2.conjugate() * P1


# #proper concatenation:
#   R = R1 * R2.conjugate()
#   print(R)
# 
# # "as-if-infinitesimal"
#   angle.inf = (w1 - w2)
# 
#   mep = w1 * e1 - w2 * e2
#   axis.inf = mep / abs(mep)
#   
  



# 
# 
# # Rotation
# u = np.array([1, 0, 0])
# R1 * u * R1.conjugate()
# 
# 
# quaternion.rotate_vectors(R=R1, v = u)
# quaternion.as_rotation_matrix(R1)
# 
# # Concatenation ("Summation of rotations")
# R1 * R2 
# # np.dot(R1, R2)
# 
# 
# w = 2 * math.acos( (p_sc * q_sc) - np.dot(p_vec, q_vec)) # angle of rotation R1 * R2
# 1 / math.sin(w/2) * (p_sc * q_vec + q_sc * p_vec + np.cross(p_vec, q_vec)) # axis of rotation R1 * R2
# 
# w == 2 * math.acos( (R2.w * R1.w) - np.dot(quaternion.as_vector_part(R2), quaternion.as_vector_part(R1)))
# 
# 
# q1_vec
# quaternion.as_vector_part(R1)
# 
# R1.w
# R1.x
# R1.y
# R1.z
# 
# quaternion.as_float_array(R1)
# quaternion.as_float_array(R1)
