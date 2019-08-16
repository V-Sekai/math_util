extends EditorPlugin
tool

func get_name() -> String: 
	return "MathUtil"

func _enter_tree() -> void:
	add_autoload_singleton("GodotMathExtension", "res://addons/math_util/math_funcs.gd")

func _exit_tree() -> void:
	remove_autoload_singleton("GodotMathExtension")
