@tool
extends Node

# Collection of smoothing methods to be used with lerp

static func spherical_to_local_position(p_theta : float, p_phi : float) -> Vector3:
	var res : Vector3 = Vector3()

	var sin_theta : float = sin(p_theta)
	var cos_theta : float = cos(p_theta)
	var sin_phi : float = sin(p_phi)
	var cos_phi : float = cos(p_phi)

	res.z = sin_theta * cos_phi
	res.y = cos_theta
	res.x = sin_theta * sin_phi
	
	return res

static func quat_from_radians(p_radians : Vector3) -> Quaternion:
	var pitch_radians : float = p_radians.x * 0.5
	var yaw_radians : float = p_radians.y * 0.5
	var roll_radians : float = p_radians.z * 0.5

	var sin_pitch : float = sin(pitch_radians)
	var cos_pitch : float = cos(pitch_radians)

	var sin_yaw : float = sin(yaw_radians)
	var cos_yaw : float = cos(yaw_radians)

	var sin_roll : float = sin(roll_radians)
	var cos_roll : float = cos(roll_radians)

	return Quaternion(sin_yaw * cos_pitch * sin_roll + cos_yaw * sin_pitch * cos_roll, sin_yaw * cos_pitch * cos_roll - cos_yaw * sin_pitch * sin_roll, cos_yaw * cos_pitch * sin_roll - sin_yaw * sin_pitch * cos_roll, cos_yaw * cos_pitch * cos_roll + sin_yaw * sin_pitch * sin_roll)

static func ease_in(t : float) -> float:
	return sin(t * PI * 0.5);
	
static func ease_out(t : float) -> float:
	return cos(t * PI * 0.5);
	
static func exponetial(t : float) -> float:
	return t * t
	
static func smooth_step(t : float) -> float:
	return t * t * (3.0 - 2.0 * t)
	
static func smoother_step(t : float) -> float:
	return t * t * t * (t * (6.0 * t - 15.0) + 10.0)
		
static func smooth_damp_scaler(current : float, target : float, current_velocity : float, smooth_time : float, max_speed : float, delta : float) -> Dictionary:
	smooth_time = max(0.0001, smooth_time)
	var value_a : float = 2.0 / smooth_time
	var value_b : float = value_a * delta
	var value_c : float = 1.0 / (1.0 + value_b + 0.48 * value_b * value_b + 0.235 * value_b * value_b * value_b)
	
	var scaler_a : float = current - target;
	var scaler_b : float = target;
	var max_length : float = max_speed * smooth_time
	
	scaler_a = clamp(scaler_a, -max_length, max_length);
	
	target = current - scaler_a;
	var scaler_c : float = (current_velocity + value_a * scaler_a) * delta;
	current_velocity = (current_velocity - value_a * scaler_c) * value_c;
	var scaler_d : float = target + (scaler_a + scaler_c) * value_c;
	if ((scaler_b - current) > 0.0 == (scaler_d > scaler_b)):
		scaler_d = scaler_b
		current_velocity = (scaler_d - scaler_b) / delta;
	
	return {"interpolation":scaler_d, "velocity":current_velocity}
		
static func smooth_damp_vector(current : Vector3, target : Vector3, current_velocity : Vector3, smooth_time : float, max_speed : float, delta : float) -> Dictionary:
	smooth_time = max(0.0001, smooth_time)
	var value_a : float = 2.0 / smooth_time
	var value_b : float = value_a * delta
	var value_c : float = 1.0 / (1.0 + value_b + 0.48 * value_b * value_b + 0.235 * value_b * value_b * value_b)
	
	var vector_a : Vector3 = current - target
	var vector_b : Vector3 = target
	var max_length = max_speed * smooth_time
	
	if (vector_a.length_squared() > max_length * max_length):
		vector_a =  vector_a.normalized() * max_length
	
	target = current - vector_a
	var vector_c : Vector3 = (current_velocity + (vector_a * value_a)) * delta
	current_velocity = (current_velocity - (vector_c * value_a)) * value_c
	var vector_d : Vector3 = target + (vector_a + vector_c) * value_c
	
	if ((vector_b - current).dot(vector_d - vector_d) > 0.0):
		vector_d = vector_b
		current_velocity = (vector_d - vector_b) / delta
		
	return {"interpolation":vector_d, "velocity":current_velocity}
	
static func camera_get_position_distance(p_camera : Camera3D, p_pos : Vector3) -> float:
	var t : Transform3D = p_camera.get_global_transform();
	var axis : Vector3 =  -Vector3(t.basis.z.x, t.basis.z.y, t.basis.z.z)
	var eyedir : Vector3 = axis.normalized()
	return eyedir.dot(p_pos) - (eyedir.dot(t.origin))
	
static func get_2d_position_from_3d_position_with_screen_limits(camera : Camera3D, position_3d : Vector3, screen_size : Vector2, screen_center : Vector2, screen_mins : Vector2, screen_max : Vector2) -> Vector2:
	var is_behind : bool = camera_get_position_distance(camera, position_3d) < 0
	var screen_pos : Vector2 = camera.unproject_position(position_3d)
	
	var screen_bounds_min : Vector2 = screen_center - screen_mins
	var screen_bounds_max : Vector2 = screen_center - (screen_size - screen_max)
	
	if !is_behind and \
		(screen_pos.x > (screen_mins.x) and screen_pos.x < (screen_max.x)) and \
		(screen_pos.y > (screen_mins.y) and screen_pos.y < (screen_max.y)):
		pass
	else:
		var rotation : int = 270
		if is_behind:
			rotation = 90
		else:
			rotation = 270
		
		screen_pos = screen_pos - screen_center
		
		var angle : float = atan2(screen_pos.y, screen_pos.x)
		angle = angle - rotation * ((PI * 2) / 360)
		
		var angle_cos : float = cos(angle)
		var angle_sin : float = sin(angle)
		
		var m : float = angle_cos / angle_sin
			
		if(angle_cos > 0):
			screen_pos = Vector2(screen_bounds_min.y / m, screen_bounds_min.y)
		else:
			screen_pos = Vector2(-screen_bounds_max.y / m, -screen_bounds_max.y)
		
		if(screen_pos.x > screen_bounds_max.x):
			screen_pos = Vector2(screen_bounds_max.x, screen_bounds_max.x*m)
		elif(screen_pos.x < -screen_bounds_min.x):
			screen_pos = Vector2(-screen_bounds_min.x, -screen_bounds_min.x*m)

		screen_pos = screen_pos + screen_center
		screen_pos.y = screen_size.y - screen_pos.y
	
	return screen_pos
	
static func get_2d_position_from_3d_position(camera : Camera3D, position_3d : Vector3) -> Vector2:
	return camera.unproject_position(position_3d)
	
static func clamp_angle(val : float, ang_min : float, ang_max : float) -> float:
	if (val < -360):
		val += 360
	if (val > 360):
		val -= 360
	
	return clamp(val, ang_min, ang_max)
	
static func adjust_facing(p_facing : Vector3, p_target : Vector3, p_step : float, p_adjust_rate : float, current_gn : Vector3):
	var n : Vector3 = p_target # normal
	var t : Vector3 = n.cross(current_gn).normalized()
	
	var x : float = n.dot(p_facing)
	var y : float = t.dot(p_facing)
	
	var ang : float = atan2(y,x)
	
	if (abs(ang) < 0.001): # too small
		return p_facing
	
	var s : float = sign(ang)
	ang = ang * s
	var turn : float = ang * p_adjust_rate * p_step
	var a : float
	if (ang<turn):
		a=ang
	else:
		a=turn
	ang = (ang - a) * s
	
	return ((n * cos(ang)) + (t * sin(ang))) * p_facing.length()
	
static func rotate_around(p_transform : Transform3D, p_point : Vector3, p_axis : Vector3, p_angle : float) -> Transform3D:
	var vector : Vector3 = p_point + (Quaternion(p_axis, p_angle) * (p_transform.origin - p_point))
	p_transform.origin = vector
	
	return p_transform.rotated(p_axis, p_angle * 0.0174532924)


static func base_log(a : float, new_base : float) -> float:
	return log(a) / log(new_base)
	
static func transform_directon_vector(p_direction : Vector3, p_basis : Basis) -> Vector3:
	return Vector3(((p_basis.x.x * p_direction.x) + (p_basis.y.x * p_direction.y) + (p_basis.z.x * p_direction.z)), ((p_basis.x.y * p_direction.x) + (p_basis.y.y * p_direction.y) + (p_basis.z.y * p_direction.z)),((p_basis.x.z * p_direction.x) + (p_basis.y.z * p_direction.y) + (p_basis.z.z * p_direction.z)))

static func get_interpolated_transform(p_current_transform : Transform3D, p_target_transform : Transform3D, p_origin_interpolation_factor : float, p_rotation_interpolation_factor : float, p_delta : float):
	var current_origin : Vector3 = p_current_transform.origin
	var current_rotation : Basis = p_current_transform.basis
	
	var target_origin : Vector3 = p_target_transform.origin
	var target_rotation : Basis = p_target_transform.basis
	
	if p_origin_interpolation_factor > 0.0:
		current_origin = current_origin.lerp(target_origin, p_origin_interpolation_factor * p_delta)
	else:
		current_origin = target_origin
		
	if p_rotation_interpolation_factor > 0.0:
		current_rotation = current_rotation.slerp(target_rotation, p_rotation_interpolation_factor * p_delta)
	else:
		current_rotation = target_rotation
	
	return Transform3D(current_rotation, current_origin)

static func vec2rad2deg(p_vector2 : Vector2) -> Vector2:
	return Vector2(rad2deg(p_vector2.x), rad2deg(p_vector2.y))

static func vec2deg2rad(p_vector2 : Vector2) -> Vector2:
	return Vector2(deg2rad(p_vector2.x), deg2rad(p_vector2.y))

static func vec3rad2deg(p_vector3 : Vector3) -> Vector3:
	return Vector3(rad2deg(p_vector3.x), rad2deg(p_vector3.y), rad2deg(p_vector3.z))

static func vec3deg2rad(p_vector3 : Vector3) -> Vector3:
	return Vector3(deg2rad(p_vector3.x), deg2rad(p_vector3.y), deg2rad(p_vector3.z))

static func sanitise_float(p_float : float) -> float:
	var return_float : float = p_float
	if is_nan(return_float) or is_inf(return_float):
		return_float = 0.0
		
	return return_float

static func sanitise_vec2(p_vec2 : Vector2) -> Vector2:
	var return_vec2 : Vector2 = p_vec2
	if is_nan(return_vec2.x) or is_inf(return_vec2.x):
		return_vec2.x = 0.0
	if is_nan(return_vec2.y) or is_inf(return_vec2.y):
		return_vec2.y = 0.0
		
	return return_vec2

static func sanitise_vec3(p_vec3 : Vector3) -> Vector3:
	var return_vec3 : Vector3 = p_vec3
	if is_nan(return_vec3.x) or is_inf(return_vec3.x):
		return_vec3.x = 0.0
	if is_nan(return_vec3.y) or is_inf(return_vec3.y):
		return_vec3.y = 0.0
	if is_nan(return_vec3.z) or is_inf(return_vec3.z):
		return_vec3.z = 0.0
		
	return return_vec3
	
static func sanitise_quat(p_quat : Quaternion) -> Quaternion:
	var return_quat : Quaternion = p_quat.normalized()
	if is_nan(return_quat.x) or is_inf(return_quat.x):
		return_quat.x = 0.0
	if is_nan(return_quat.y) or is_inf(return_quat.y):
		return_quat.y = 0.0
	if is_nan(return_quat.z) or is_inf(return_quat.z):
		return_quat.z = 0.0
	if is_nan(return_quat.w) or is_inf(return_quat.w):
		return_quat.w = 1.0
		
	return return_quat
	
static func shortest_angle_distance(p_from: float, p_to: float) -> float:
	var difference: float = fmod(p_to - p_from, PI * 2.0)
	return fmod(2.0 * difference, PI * 2.0) - difference

static func minimum_velocity(p_velocity: Vector3, p_min: float) -> Vector3:
	if abs(p_velocity.x) <= p_min:
		p_velocity.x = 0.0
		
	if abs(p_velocity.y) <= p_min:
		p_velocity.y = 0.0
		
	if abs(p_velocity.z) <= p_min:
		p_velocity.z = 0.0
		
	return p_velocity

static func degrade_velocity(p_velocity: Vector3, p_degrade: float, p_min: float) -> Vector3:
	if abs(p_velocity.x) + p_min > p_degrade:
		p_velocity.x -= p_degrade * sign(p_velocity.x)
	else:
		p_velocity.x = 0.0
		
	if abs(p_velocity.y) + p_min > p_degrade:
		p_velocity.y -= p_degrade * sign(p_velocity.y)
	else:
		p_velocity.y = 0.0
		
	if abs(p_velocity.z) + p_min > p_degrade:
		p_velocity.z -= p_degrade * sign(p_velocity.z)
	else:
		p_velocity.z = 0.0
		
	return p_velocity
