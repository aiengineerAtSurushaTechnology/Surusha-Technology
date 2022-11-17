#ifndef SCENE_H
#define SCENE_H

#include <QGraphicsScene>

class Scene : public QGraphicsScene
{
    Q_OBJECT
public:
    enum ToolType {
        Cursor,
        Pen,
        Rect,
        Ellipse,
        Star,
        QtQuick,
        Eraser
    };
    explicit Scene(QObject *parent = nullptr);

    ToolType getTool() const;
    void setTool(ToolType newTool);

signals:

public slots:

    // QGraphicsScene interface
protected:
    void dragMoveEvent(QGraphicsSceneDragDropEvent *event) override;
    void dropEvent(QGraphicsSceneDragDropEvent *event) override;



    // QGraphicsScene interface
protected:
    void mousePressEvent(QGraphicsSceneMouseEvent *event) override;
    void mouseMoveEvent(QGraphicsSceneMouseEvent *event) override;
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event) override;

private:

    void drawLineTo(const QPointF &endPoint);
    void drawEraserAt(const QPointF & endPoint);
    void eraseStrokeUnder(QGraphicsEllipseItem * item);

    ToolType tool;
    bool drawing;
    QGraphicsItemGroup * lineGroup;
    QPointF startingPoint;
    QPointF lastPenPoint;
    QGraphicsEllipseItem * lastEraserCircle;

    QGraphicsLineItem * horGuideLine;
    QGraphicsLineItem * verGuideLine;
};

#endif // SCENE_H
