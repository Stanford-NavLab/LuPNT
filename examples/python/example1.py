from pylupnt import hello, VectorDouble, multiply_by_scalar


def main():
    print(hello("Jeremy"))

    vd = VectorDouble([0.0, 1.0, 2.0, 3.0])
    vd2 = multiply_by_scalar(vd, 2.0)
    print(vd)
    print("Multiplied by 2.0 is:")
    print(vd2)


if __name__ == "__main__":
    main()
