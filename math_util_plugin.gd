extends EditorPlugin
tool

func get_name(): 
	return "MathUtil"

func _enter_tree():
	add_autoload_singleton("GodotMathExtension", "res://addons/math_util/math_funcs.gd")

func _exit_tree():
	remove_autoload_singleton("GodotMathExtension")
