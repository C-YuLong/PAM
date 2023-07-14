set files [exec sh -c {find -name "*.v" -type f}]
foreach file $files {
    analyze -format verilog -library WORK $file
}
elaborate $name -architecture "verilog" -library WORK