# Go to the root of the project
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# Configure and build
cmake -Sexamples/cpp -Bbuild_examples -DCMAKE_BUILD_TYPE=Release
cmake --build build_examples -j4 --target all_examples
# Run
cd ./build_examples/examples
for example in *; do
    file_path=$(find ../../examples/cpp -name "$example.cc")
    if [ -n "$file_path" ]; then
        if ! grep -qE "show()|draw()" "$file_path"; then
            echo -e "\n\n*********** Running example: $example ***********\n"
            ./"$example"
            exit_status=$?
            if [ $exit_status -ne 0 ]; then
                echo -e "\n\n********** FAILED **********"
                exit 1
            else
                echo -e "\n\n********** SUCCESS **********"
            fi
        fi
    fi
done
cd ..