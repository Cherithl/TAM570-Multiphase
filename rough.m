clc
close all

[X_Mat] = stitch(X, Nelx,Nely, N);
[Y_Mat] = stitch(Y, Nelx,Nely, N);
[U_Mat] = stitch(U, Nelx,Nely, N);
[V_Mat] = stitch(V, Nelx,Nely, N);
[T_Mat] = stitch(T, Nelx,Nely, N);


figure
    contourf(X_Mat, Y_Mat, T_Mat, 80, "Linecolor","None")
    colorbar
    colormap("hot")
    hold on
    quiver(X_Mat,Y_Mat,U_Mat,V_Mat,2,'LineWidth',2,'Color', 'k'); hold on;
    % quiver(X_Mat,Y_Mat,U_Mat,V_Mat,'off','Color', 'k'); hold on;
    axis equal

figure
    mesh(X_Mat, Y_Mat, zeros(size(X_Mat)),'LineWidth',2, 'EdgeColor','k')
    axis equal
    view(2)