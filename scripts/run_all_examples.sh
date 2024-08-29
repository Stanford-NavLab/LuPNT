# Go to the root of the project
ROOT_DIR=$(dirname "$(realpath "$0")")/..
cd $ROOT_DIR
# Configure and build
cmake -Sexamples/cpp -Bbuild_examples -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build_examples -j12 --target all_examples
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
                cd $ROOT_DIR
                exit 1
            else
                echo -e "\n\n********** SUCCESS **********"
            fi
        fi
    fi
done
cd $ROOT_DIR
